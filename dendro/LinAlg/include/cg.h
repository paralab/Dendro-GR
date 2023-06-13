//
// Created by milinda on 12/6/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This is modified to perform matrix-free iterative methods for solving linear systems that
 * generated from the Galerkin-Discretization (FEM).
*/
//

#ifndef SFCSORTBENCH_CG_H
#define SFCSORTBENCH_CG_H

#include "mathUtils.h"
#include "mesh.h"
#include "operators.h"

// Code taken from : http://math.nist.gov/iml++/
//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//
// Note: The above method is modified for simple, CG without preconditioning.
//
//*****************************************************************
namespace linalg
{

    template <typename Real >
    int CG(ot::Mesh* mesh, Real *x,Real *b, int max_iter, Real& tol)
    {

        Real resid,alpha,beta,rho,rho_1;
        int status=1; // 0 indicates it has solved the system within the specified max_iter, 1 otherwise.

        if(mesh->isActive())
        {

            int activeRank=mesh->getMPIRank();
            int activeNpes=mesh->getMPICommSize();
            MPI_Comm activeComm=mesh->getMPICommunicator();

            const unsigned int dof=mesh->getDegOfFreedom();
            const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
            const unsigned int local_dof=nodeLocalEnd-nodeLocalBegin;


            // need to release the memory before exit.
            Real * p=new Real[dof];
            Real * z=new Real[dof];
            Real * q=new Real[dof];

            Real * Ax=new Real[dof];
            Real * Ap=new Real[dof];
            Real * r0=new Real[dof];
            Real * r1=new Real[dof];


            Real normb = normLInfty(b+nodeLocalBegin,local_dof,activeComm);
            par::Mpi_Bcast(&normb,1,0,activeComm);
            fem::operators::poisson::matvec(mesh,x,Ax);

            char fPrefix[256];
            sprintf(fPrefix,"%s_%d","cg",0);
            const char * varNames[]={"U"};
            const double * var[]={Ax};
            io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);

            for(unsigned int i=nodeLocalBegin;i<nodeLocalEnd;i++)
            {
                r0[i]=b[i]-Ax[i];
                p[i]=r0[i];
            }


            if (normb == 0.0)
                normb = 1;

            Real normr=normLInfty(r0+nodeLocalBegin,local_dof,activeComm);
            par::Mpi_Bcast(&normr,1,0,activeComm);
            if(!activeRank) std::cout<<"initial residual : "<<(normr/normb)<<std::endl;

            if ((resid = normr / normb) <= tol) {
                tol = resid;
                max_iter = 0;

                delete [] p;
                delete [] z;
                delete [] q;
                delete [] Ax;
                delete [] Ap;
                delete [] r0;
                delete [] r1;

                status=0;
            }

            if(status!=0)
            {

                for(unsigned int i=1;i<=max_iter;i++)
                {
                    fem::operators::poisson::matvec(mesh,p,Ap, false);

                    alpha=(dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)/dot(p+nodeLocalBegin,Ap+nodeLocalBegin,local_dof,activeComm));
                    par::Mpi_Bcast(&alpha,1,0,activeComm);

                    //if(!activeRank) std::cout<<"<r_0,r_0> : "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" <p,Ap>: "<<dot(p+nodeLocalBegin,Ap+nodeLocalBegin,local_dof,activeComm)<<" alpha "<<alpha<<std::endl;

                    /*Real normAx=normL2(Ax+nodeLocalBegin,local_dof,activeComm);
                    Real normx=normL2(x+nodeLocalBegin,local_dof,activeComm);
                    if(!activeRank)std::cout<<"l2(Ax): "<<normAx<<" l2(x): "<<normx<<std::endl;*/

                    //if(!activeRank) std::cout<<"rank: " <<activeRank<<" alpha: "<<alpha<<std::endl;
                    //std::cout<<"alpha: "<<alpha<<std::endl;

                    for(unsigned int e=nodeLocalBegin;e<nodeLocalEnd;e++)
                    {
                        x[e]+=alpha*p[e];
                        r1[e]=r0[e]-alpha*Ap[e];
                    }

                    normr=normLInfty(r1+nodeLocalBegin,local_dof,activeComm);
                    par::Mpi_Bcast(&normr,1,0,activeComm);

                    if((!activeRank) && (i%10==0)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;

                    if ((resid = normr / normb) <= tol) {

                        if((!activeRank)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;
                        tol = resid;
                        max_iter = i;

                        delete [] p;
                        delete [] z;
                        delete [] q;
                        delete [] Ax;
                        delete [] Ap;
                        delete [] r0;
                        delete [] r1;

                        status=0;
                        break;
                    }

                    beta=(dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)/dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm));
                    par::Mpi_Bcast(&beta,1,0,activeComm);

                    //if(!activeRank) std::cout<<"<r_1,r_1> : "<<dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)<<" <r_0,r_0>: "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" beta "<<beta<<std::endl;



                    for(unsigned int e=nodeLocalBegin;e<nodeLocalEnd;e++)
                    {
                        p[e]=r1[e]+beta*p[e];
                        r0[e]=r1[e];
                    }

                    fem::operators::poisson::applyDirichletBoundary(mesh,p,0.0);
                    fem::operators::poisson::applyDirichletBoundary(mesh,r0,0.0);


                    char fPrefix[256];
                    sprintf(fPrefix,"%s_%d","cg",i);
                    const char * varNames[]={"U"};
                    const double * var[]={x};
                    io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);


                }

                if(status!=0)
                {
                    tol = resid;
                    delete [] p;
                    delete [] z;
                    delete [] q;
                    delete [] Ax;
                    delete [] r0;
                    delete [] r1;
                    status=1;

                }



            }


        }

        // bcast act as a barrier for active and inactive meshes.
        par::Mpi_Bcast(&tol,1,0,mesh->getMPIGlobalCommunicator());
        return status;

    }


}



#endif //SFCSORTBENCH_CG_H
