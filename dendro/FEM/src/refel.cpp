//
// Created by milinda on 12/26/16.
//

/**
 * @author Milinda Fernando
 * @author Hari Sundar
 * @breif Contains data structures to store the reference element information.
 *
 * @refference: Based of HOMG code written in matlab.
 * */

#include "refel.h"


RefElement::RefElement()// default constructor
{

}

RefElement::RefElement(unsigned int dim, unsigned int order)
{

    /*
     * Reference element domain is from  (-1,1)
     * */
    if(dim==1)
        m_uiVol=2;
    else if(dim==2)
        m_uiVol=4;
    else if(dim==3)
        m_uiVol=8;

    m_uiDimension=dim;
    m_uiOrder=order;

    m_uiNp  = (m_uiOrder + 1) * (m_uiOrder + 1) * (m_uiOrder + 1);
    m_uiNfp = (m_uiOrder + 1) * (m_uiOrder + 1);
    m_uiNrp = (m_uiOrder + 1);

    u.resize(m_uiNrp);
    r.resize(m_uiNrp);
    g.resize(m_uiNrp);
    u_0.resize(m_uiNrp);
    u_1.resize(m_uiNrp);

    Vu.resize(m_uiNrp*m_uiNrp);
    gradVu.resize(m_uiNrp*m_uiNrp);

    w.resize(m_uiNrp);
    wgll.resize(m_uiNrp);

    Vr.resize(m_uiNrp*m_uiNrp);
    gradVr.resize(m_uiNrp*m_uiNrp);

    Vg.resize(m_uiNrp*m_uiNrp);
    gradVg.resize(m_uiNrp*m_uiNrp);

    Vu_0.resize(m_uiNrp*m_uiNrp);
    Vu_1.resize(m_uiNrp*m_uiNrp);

    ip_1D_0.resize(m_uiNrp*m_uiNrp);
    ip_1D_1.resize(m_uiNrp*m_uiNrp);

    quad_1D.resize(m_uiNrp*m_uiNrp);
    quadT_1D.resize(m_uiNrp*m_uiNrp);

    Dr.resize(m_uiNrp*m_uiNrp);
    Dg.resize(m_uiNrp*m_uiNrp);
    DgT.resize(m_uiNrp*m_uiNrp);

    im_vec1.resize(m_uiNp); // resize intermidiate values for number of points in 1D case. 
    im_vec2.resize(m_uiNp); // resize intermidiate values for number of points in 1D case. 

    ipT_1D_0.resize(m_uiNrp*m_uiNrp);
    ipT_1D_1.resize(m_uiNrp*m_uiNrp);
    

    #ifdef WITH_BLAS_LAPACK

        unsigned int info;

        double x_min=-1.0;
        double x_max= 1.0;

        for(unsigned int k=0;k<m_uiNrp;k++){
            u[k]=x_min+(x_max-x_min)*k/m_uiOrder;
        }

        //std::cout<<"gauss-labatto quadrature begin "<<std::endl;

        double* w1=new double[m_uiNrp];
        double* wgll1=new double[m_uiNrp];

        const double alpha =0.0;
        const double beta =0.0;

        basis::jacobiglq(alpha,beta,m_uiOrder,&(*(r.begin())),wgll1);
        basis::jacobigq(alpha,beta,m_uiOrder,&(*(g.begin())),w1);
        //std::cout<<"gauss-labatto quadrature end "<<std::endl;

        for(unsigned int k=0;k<m_uiNrp;k++)
        {
            // parent to child points.
            u_0[k]=0.5*(u[k]-1);
            u_1[k]=0.5*(u[k]+1);
        }

        /*std::cout<<" r: ";printArray_1D((&(*(r.begin()))),m_uiNrp);
        std::cout<<" u: ";printArray_1D((&(*(u.begin()))),m_uiNrp);
        std::cout<<" wgll: ";printArray_1D(wgll1,m_uiNrp);
        std::cout<<" wg: ";printArray_1D(w1,m_uiNrp);
        std::cout<<" g: ";printArray_1D((&(*(g.begin()))),m_uiNrp);

        std::cout<<" r0: ";printArray_1D((&(*(u_0.begin()))),m_uiNrp);
        std::cout<<" r1: ";printArray_1D((&(*(u_1.begin()))),m_uiNrp);*/


        for(unsigned int i=0;i<m_uiNrp;i++) {


            basis::jacobip(alpha,beta, i, (&(*(u_0.begin()))), &(*(Vu_0.begin() + i*m_uiNrp )), m_uiNrp);
            basis::jacobip(alpha,beta, i, (&(*(u_1.begin()))), &(*(Vu_1.begin() + i*m_uiNrp )), m_uiNrp);

            basis::jacobip(alpha,beta, i, (&(*(u.begin()))), &(*(Vu.begin() + i*m_uiNrp)), m_uiNrp);
            basis::gradjacobip(alpha,beta, i, (&(*(u.begin()))), &(*(gradVu.begin() + i*m_uiNrp)), m_uiNrp);

            basis::jacobip(alpha,beta, i, (&(*(r.begin()))), &(*(Vr.begin() + i*m_uiNrp)), m_uiNrp);
            basis::gradjacobip(alpha,beta, i, (&(*(r.begin()))), &(*(gradVr.begin() + i*m_uiNrp)), m_uiNrp);

            basis::jacobip(alpha,beta, i, (&(*(g.begin()))), &(*(Vg.begin() + i*m_uiNrp)), m_uiNrp);
            basis::gradjacobip(alpha,beta, i, (&(*(g.begin()))), &(*(gradVg.begin() + i*m_uiNrp)), m_uiNrp);

            //std::cout<<i<<" r eval : ";printArray_1D(&(*(Vr.begin() + i*m_uiNrp)),m_uiNrp);
            //std::cout<<i<<" r0 eval : ";printArray_1D(&(*(xCh_0.begin())),m_uiNrp);
            //std::cout<<i<<" r1:eval  ";printArray_1D(&(*(xCh_1.begin())),m_uiNrp);

        }

            /*std::cout<<"Vr0: "<<std::endl;
            printArray_2D(&(*(Vu_0.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"Vr1: "<<std::endl;
            printArray_2D(&(*(Vu_1.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"Vu: "<<std::endl;
            printArray_2D(&(*(Vu.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"Vr: "<<std::endl;
            printArray_2D(&(*(Vr.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"Vg: "<<std::endl;
            printArray_2D(&(*(Vg.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"gradVu: "<<std::endl;
            printArray_2D(&(*(gradVu.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"gradVr: "<<std::endl;
            printArray_2D(&(*(gradVr.begin())),m_uiNrp,m_uiNrp);

            std::cout<<"gradVg: "<<std::endl;
            printArray_2D(&(*(gradVg.begin())),m_uiNrp,m_uiNrp);*/


            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu.begin())),m_uiNrp,&(*(gradVr.begin())),&(*(Dr.begin())),m_uiNrp,info);
            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu.begin())),m_uiNrp,&(*(gradVg.begin())),&(*(Dg.begin())),m_uiNrp,info);


            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu.begin())),m_uiNrp,&(*(Vu_0.begin())),&(*(ip_1D_0.begin())),m_uiNrp,info);
            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu.begin())),m_uiNrp,&(*(Vu_1.begin())),&(*(ip_1D_1.begin())),m_uiNrp,info);

            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu.begin())),m_uiNrp,&(*(Vg.begin())),&(*(quad_1D.begin())),m_uiNrp,info);


            // transpose operators
            /*lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu_0.begin())),m_uiNrp,&(*(Vu.begin())),&(*(ipT_1D_0.begin())),m_uiNrp,info);
            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vu_1.begin())),m_uiNrp,&(*(Vu.begin())),&(*(ipT_1D_1.begin())),m_uiNrp,info);

            lapack::lapack_DGESV(m_uiNrp,m_uiNrp,&(*(Vg.begin())),m_uiNrp,&(*(Vu.begin())),&(*(quadT_1D.begin())),m_uiNrp,info);*/



            for(unsigned int i=0;i<m_uiNrp;i++) {
                w[i]=w1[i];
                wgll[i]=wgll1[i];

                for(unsigned int j=0;j<m_uiNrp;j++)
                {
                    ipT_1D_0[i*m_uiNrp+j]=ip_1D_0[j*m_uiNrp+i];
                    ipT_1D_1[i*m_uiNrp+j]=ip_1D_1[j*m_uiNrp+i];
                    quadT_1D[i*m_uiNrp+j]=quad_1D[j*m_uiNrp+i];
                    DgT[i*m_uiNrp+j]=Dg[j*m_uiNrp+i];


                }

            }


            delete [] w1;
            delete [] wgll1;

        /*std::cout<<"ip_1D_0: "<<std::endl;
        printArray_2D(&(*(ip_1D_0.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"ip_1D_1: "<<std::endl;
        printArray_2D(&(*(ip_1D_1.begin())),m_uiNrp,m_uiNrp);


        std::cout<<"ipT_1D_0: "<<std::endl;
        printArray_2D(&(*(ipT_1D_0.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"ipT_1D_1: "<<std::endl;
        printArray_2D(&(*(ipT_1D_1.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"Dr: "<<std::endl;
        printArray_2D(&(*(Dr.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"Dg: "<<std::endl;
        printArray_2D(&(*(Dg.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"Q: "<<std::endl;
        printArray_2D(&(*(quad_1D.begin())),m_uiNrp,m_uiNrp);

        std::cout<<"QT: "<<std::endl;
        printArray_2D(&(*(quadT_1D.begin())),m_uiNrp,m_uiNrp);*/
        //std::cout<<" wg: ";printArray_1D((&(*(w.begin()))),m_uiNrp);

        //printArray_2D(&(*(Fr.begin())),m_uiNrp,m_uiNrp);
    
    #else


        if(m_uiDimension==3 && m_uiOrder==1)
        {
            for(unsigned int i=0;i<m_uiNrp;i++ )
            {
            for(unsigned int j=0;j<m_uiNrp;j++)
            {
                ip_1D_0[i*m_uiNrp+j]=IP_1D_Order_1_0[i*m_uiNrp+j];
                ipT_1D_0[i*m_uiNrp+j]=IP_1D_Order_1_0[j*m_uiNrp+i];

                ip_1D_1[i*m_uiNrp+j]=IP_1D_Order_1_1[i*m_uiNrp+j];
                ipT_1D_1[i*m_uiNrp+j]=IP_1D_Order_1_1[j*m_uiNrp+i];
            }

            }

        }else if (m_uiDimension==3 && m_uiOrder==2)
        {
            for(unsigned int i=0;i<m_uiNrp;i++ )
            {
                for(unsigned int j=0;j<m_uiNrp;j++)
                {
                    ip_1D_0[i*m_uiNrp+j]=IP_1D_Order_2_0[i*m_uiNrp+j];
                    ipT_1D_0[i*m_uiNrp+j]=IP_1D_Order_2_0[j*m_uiNrp+i];

                    ip_1D_1[i*m_uiNrp+j]=IP_1D_Order_2_1[i*m_uiNrp+j];
                    ipT_1D_1[i*m_uiNrp+j]=IP_1D_Order_2_1[j*m_uiNrp+i];

                }

            }

        }else if (m_uiDimension==3 && m_uiOrder==4)
        {
            for(unsigned int i=0;i<m_uiNrp;i++ )
            {
                for(unsigned int j=0;j<m_uiNrp;j++)
                {
                    ip_1D_0[i*m_uiNrp+j]=IP_1D_Order_4_0[i*m_uiNrp+j];
                    ipT_1D_0[i*m_uiNrp+j]=IP_1D_Order_4_0[j*m_uiNrp+i];

                    ip_1D_1[i*m_uiNrp+j]=IP_1D_Order_4_1[i*m_uiNrp+j];
                    ipT_1D_1[i*m_uiNrp+j]=IP_1D_Order_4_1[j*m_uiNrp+i];
                }

            }

        }else
        {
            std::cout<<"RefEl: Error invalid dimension and order specified"<<std::endl;
        }


    #endif

     // computes the p2c and c2p operators.
    p2c.resize((2*m_uiOrder+1)*m_uiNrp,0.0);
    c2p.resize((2*m_uiOrder+1)*m_uiNrp,0.0);

    // child 0  and child 1
    for(unsigned int i=0; i < m_uiNrp ; i++)
    {
        for(unsigned int j=0; j < m_uiNrp ; j++)
        {
            p2c[(0*m_uiOrder + i )*m_uiNrp + j ] = ipT_1D_0[i*m_uiNrp+j];
            p2c[(1*m_uiOrder + i )*m_uiNrp + j ] = ipT_1D_1[i*m_uiNrp+j];
        }

    }

    for(unsigned int i=0; i < (2*m_uiOrder+1) ; i++)
        for(unsigned int j=0; j < m_uiNrp ; j++)
            c2p[j*(2*m_uiOrder+1) + i ] = p2c[i*(m_uiNrp) + j]; 

    //std::cout<<"ele order : "<<m_uiOrder<<" p2c"<<std::endl;
    //printArray_2D(p2c.data(),(2*m_uiOrder+1),m_uiNrp);
    //printArray_2D(c2p.data(),m_uiNrp,(2*m_uiOrder+1));

    const unsigned int npb2 = 2*m_uiOrder+1;
    c2p_2d.resize((m_uiNrp*m_uiNrp) * (npb2*npb2));
    c2p_3d.resize((m_uiNrp*m_uiNrp*m_uiNrp) * (npb2*npb2*npb2) );

    kron(c2p.data(),c2p.data(), c2p_2d.data(), m_uiNrp, npb2, m_uiNrp, npb2);
    kron(c2p.data(),c2p_2d.data(), c2p_3d.data(), m_uiNrp, npb2, m_uiNrp*m_uiNrp, npb2*npb2 );

    // std::cout<<"c2p-3d: \n";
    // printArray_2D(c2p_3d.data(),m_uiNp,(npb2*npb2*npb2));
    // std::cout<<std::endl;


    
}

RefElement::~RefElement() {

    im_vec1.clear();
    im_vec2.clear();

}


void RefElement::generateHeaderFile(char * fName)
{


    #ifdef WITH_BLAS_LAPACK
        std::ofstream myfile (fName);
        if (myfile.is_open())
        {

        myfile<<" /**\n * @author Milinda Fernando\n * @breif This file contains the precomputed interpolation matrices for a given order, computed based on Hari's HOMG refele.m code. \n * @Note: For wavelet based finite differencing we use uniform points instead of gll points.\n */";
        myfile << "// This is a machine generated code to initialize interpolation matrices for reference element specified. \n";

        myfile<<"static double IP_1D_Order_"<<m_uiOrder<<"_0 [] ={";
        for(unsigned int i=0;i<m_uiNrp;i++)
            for(unsigned int j=0;j<m_uiNrp;j++)
                (i!=(m_uiNrp-1) || j!=(m_uiNrp-1)) ? myfile<<std::setprecision(16)<<ip_1D_0[i*m_uiNrp+j]<<" ," : myfile<<std::setprecision(16)<<ip_1D_0[i*m_uiNrp+j];

            myfile<<" };\n";


            myfile<<"static double IP_1D_Order_"<<m_uiOrder<<"_1 [] ={";
            for(unsigned int i=0;i<m_uiNrp;i++)
                for(unsigned int j=0;j<m_uiNrp;j++)
                    (i!=(m_uiNrp-1) || j!=(m_uiNrp-1)) ? myfile<<std::setprecision(16)<<ip_1D_1[i*m_uiNrp+j]<<" ," : myfile<<std::setprecision(16)<<ip_1D_1[i*m_uiNrp+j];

            myfile<<" };\n";
        myfile.close();
        }
        else std::cout << "Unable to open file"<<std::endl;

    #else
        std::cout<<"GenerateHeader file should be run with BLAS enabled. "<<std::endl;
    #endif

}



void RefElement::I3D_Children2ParentInjection(const double * in, double* out) const
{

    const unsigned int npb2 = 2*m_uiOrder+1;
    
    // pure injection ------
    for(unsigned int k=0; k < npb2; k++)
     for(unsigned int j=0; j < npb2; j++)
      for(unsigned int i=0; i < npb2; i++)
        if(k%2==0 && j%2==0 && i%2==0)
            out[(k>>1u)*m_uiNrp*m_uiNrp + (j>>1u)*m_uiNrp + (i>>1u) ] = in[k*npb2*npb2 + j*npb2 + i];


}


void RefElement::I3D_Children2Parent(const double * in, double* out) const
{

    double * im1=(double *)&(*(im_vec1.begin()));
    double * im2=(double *)&(*(im_vec2.begin()));

    const unsigned int npb2 = 2*m_uiOrder+1;

    const unsigned int M = m_uiNp;
    const unsigned int N = (npb2*npb2*npb2);

    for(unsigned int m=0; m < M; m++)
    {
        out[m]=0;
        for(unsigned int n=0; n < N; n++)
            out[m] += c2p_3d[m * (npb2*npb2*npb2) + n] * in[n];
    
    }


}
