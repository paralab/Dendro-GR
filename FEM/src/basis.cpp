//
// Created by milinda on 1/12/17.
//

/**
 * @author Milinda Fernando
 * School of Computing, University of Utah.
 * @breif Contains implementation of Jacobi polynomials, and Jacobi Gauss-Lobatto quadrature.
 * reference: Mangll implementation.
 * */
#include <limits>
#include <vector>
#include "basis.h"

namespace basis {



    void jacobip(double alpha, double beta, unsigned int N, double *x, double *p , unsigned int np)
    {


        unsigned int  xn = np;
        double  apb, gamma0, gamma1, isqrtgamma0, isqrtgamma1;
        double  aold, anew, bnew, h1;
        double  * pl=new double [(N+1)*xn];

        assert(N>=0 && alpha >-1.0 && beta >-1.0);
        //std::cout<<" xn: "<<xn<<std::endl;
        apb = alpha + beta;
        gamma0 = pow (2.0, (apb + 1.0)) / (apb + 1.0)* std::tgamma (alpha + 1.0) * std::tgamma (beta + 1.0) / std::tgamma (apb + 1.0);
        isqrtgamma0 = 1.0 / sqrt (gamma0);

        gamma1 = (alpha + 1.0) * (beta + 1.0) / (apb + 3.0) * gamma0;
        isqrtgamma1 = 1.0 / sqrt (gamma1);

        for(unsigned int k=0;k<xn;k++) {
            pl[k] = isqrtgamma0;
            //std::cout<<"pl[ "<<k<<"]: "<<pl[k]<<std::endl;
        }

        if(N>0) {
            for (unsigned int k = 0; k < xn; k++) {
                pl[xn + k] = ((alpha + beta + 2) * (x[k]) / 2.0 + (alpha - beta) / 2.0) / sqrt(gamma1);
            }
        }

        if (N == 0) {
           memcpy(p,pl,sizeof(double)*xn);
        }
        else if (N == 1) {
            memcpy(p,pl+xn,sizeof(double)*xn);
        }
        else {
            aold = 2.0/(2.0+alpha+beta)*sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0));
            for (unsigned int i = 0; i < N - 1; i++) {
                h1 = 2.0 *(i+1 )+ apb;
                anew =2.0/(h1+2.0)*sqrt((i+2)*(i+2+alpha+beta)*(i+2+alpha)*(i+2+beta)/(h1+1.0)/(h1+3.0));
                bnew = -(alpha * alpha - beta * beta) / h1 / (h1 + 2.0);
                for (unsigned int k = 0; k < xn; k++)
                    pl[(i + 2) * xn + k] =1.0 / anew * (-aold * pl[(i) * xn + k] + (x[k] - bnew) * pl[(i + 1) * xn + k]);

                aold=anew;
            }
            memcpy(p, pl + N * xn, sizeof(double) * xn);
         }
        delete [] pl;
        return  ;

    }

    void gradjacobip(double alpha, double beta, int N,double *x, double *dp, unsigned int np)
    {

        assert(N>=0 && alpha >-1.0 && beta >-1.0);
        if (N == 0) {
            for(unsigned int k=0;k<np;k++)
                dp[k]=0;
        }
        else {
            jacobip(alpha + 1.0, beta + 1.0,N-1,x,dp,np);
            for(unsigned int k=0;k<np;k++)
                dp[k]=dp[k]*sqrt (N * (N + alpha + beta + 1.0));


        }

        return ;

    }


    void jacobigq(double alpha, double beta, int N, double *x, double *w)
    {

#ifdef WITH_BLAS_LAPACK
        //Note: size of x and w should be (N+1)
        if(N==0)
        {
            x[0]=(alpha-beta)/(alpha+beta+2);
            w[0]=2;
            return ;
        }

        double *ji=new double[(N+1)*(N+1)];
        double *jtrans=new double[(N+1)*(N+1)];
        double * h1=new double [(N+1)];

        double * wr=new double[(N+1)];
        double * wi=new double[(N+1)];
        double * vs=new double[(N+1)*(N+1)];

        for(unsigned int k=0;k<(N+1);k++)
            h1[k]=2*k+alpha+beta;

        for(unsigned int i=0;i<(N+1);i++)
            for(unsigned int j=0;j<(N+1);j++)
            {
                if(i!=j)
                    ji[i*(N+1)+j]=0;
                else
                    ji[i*(N+1)+j]=(-0.5*(alpha*alpha-beta*beta)/(h1[j]+2))/h1[j];
            }

        for(unsigned int j=0;j<N;j++)
        { // added to the diagonal shifted by 1.
            ji[j*(N+1)+j+1]=2.0/(h1[j]+2.0)*sqrt((j+1)*(j+1+alpha+beta)*(j+1+alpha)*(j+1+beta)/(h1[j]+1.0)/(h1[j]+3.0));
        }
        if((alpha+beta)<10*std::numeric_limits<double>::epsilon())
            ji[0*(N+1)+0]=0;

        for(unsigned int i=0;i<(N+1);i++)
            for(unsigned int j=0;j<(N+1);j++)
                jtrans[i*(N+1)+j]=ji[i*(N+1)+j]+ji[j*(N+1)+i];

/*        for(unsigned int i=0;i<(N+1);i++)
        {
            for(unsigned int j=0;j<(N+1);j++)
                std::cout<<jtrans[i*(N+1)+j]<<" ";
            std::cout<<std::endl;
        }*/


//      Compute quadrature by eigenvalue solve
        unsigned int info;
        //std::cout<<" lapack eigen solve begin "<<std::endl;
        lapack::lapack_DSYEV((N+1),jtrans,(N+1),wr,vs,info);
        //[V,D] = eig(J); x = diag(D);
        //std::cout<<" lapack eigen solve end "<<std::endl;
        memcpy(x,wr,sizeof(double)*(N+1));

        for(unsigned int k=0;k<(N+1);k++)
            w[k]=vs[0*(N+1)+k]*vs[0*(N+1)+k]*pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

        /*for(unsigned int k=0;k<(N+1);k++)
            std::cout<<" w["<<k<<"]: "<<w[k]<<std::endl;*/

        delete [] ji;
        delete [] jtrans;
        delete [] h1;
        delete [] wr;
        delete [] wi;
        delete [] vs;
#endif

    }




    void jacobiglq(double alpha, double beta, int N, double *x, double *w)
    {

#ifdef WITH_BLAS_LAPACK
        if(N==1)
        {
            x[0]=-1.0;
            x[1]=1.0;

            w[0]=1.0;
            w[1]=1.0;
            return ;
        }


        if(N>2) {
            double * xint=new double [N-2+1];
            double * wq=new double [N-2+1];
            //std::cout<<" gauss quad begin "<<std::endl;
            basis::jacobigq(alpha + 1, beta + 1, (N - 2),xint,wq);
            //std::cout<<"gauss quad end "<<std::endl;
            for(unsigned int k=1;k<=(N-1);k++)
                x[k]=xint[k-1];

            x[0]=-1.0;
            x[N]=1.0;

            basis::jacobip(alpha,beta,N,x,w,(N+1));
            double adgammaN = (2.0*N + alpha + beta + 1.0) / (N * (alpha + beta + N + 1.0));

            for(unsigned int k=0;k<(N+1);k++)
            {
                w[k]=w[k]*w[k];
                w[k]=adgammaN/w[k];
            }

            w[0]=w[0]*(1.0+alpha);
            w[N]=w[N]*(1.0+ beta);

            delete [] xint;
            delete [] wq;

        }
#endif
    }





}