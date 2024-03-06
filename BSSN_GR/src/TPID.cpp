/**
 * @file TPID.cpp
 * @brief Invoke the TPID solver and write the soluition to the disk. 
 * This solution will be used by the BSSNSolver, to initialize the data.  
 * @version 0.1
 * @date 2020-08-11
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include <iostream>
#include "TwoPunctures.h"
#include "mpi.h"
#include "grUtils.h"

int main(int argc, char **argv)
{   

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc,&argv);

    int npes,rank;
    MPI_Comm_size(comm, &npes);

    if(npes>1)
    {
        std::cout<<"tpid should be executed single mpi task"<<std::endl;
        MPI_Abort(comm,0);
    }


    if(argc<3)
    {
        std::cout<<"ERROR: Not Enough Arguments!" << std::endl << "Usage: "<<argv[0]<<"paramFile numthreads"<<std::endl;
        return 0;
    }

    const unsigned int nthreads =atoi(argv[2]);
    TP_OMP_THREADS=nthreads;
    
    std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    bssn::readParamFile(argv[1],comm);

    double massp, massm, massp_adm, massm_adm, E, J1, J2, J3;
    double t1 = MPI_Wtime();
    //TwoPunctures(&cctkGH, &massp, &massm, &massp_adm, &massm_adm,&E, &J1, &J2, &J3);
    TP_SOLVE_VALID=0;
    if(TP_SOLVE_VALID==0)
    {
        TPStore(&massp, &massm, &massp_adm, &massm_adm,&E, &J1, &J2, &J3,TPID::FILE_PREFIX.c_str());
        TP_SOLVE_VALID=1;
    }               
    double t2 = MPI_Wtime();
    std::cout<<"TP solver time : "<<(t2-t1)<<std::endl;

    MPI_Finalize();
    return 0;


}
