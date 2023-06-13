/**
 * @file launcher.h
 * @author Milinda Fernando
 * @brief Class to maintain launching different solvers in different MPI communicators at the same time. 
 * @version 0.1
 * @date 2021-10-12
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once
#include<iostream>
#include<vector>
#include "dendro.h"
#include "mpi.h"


namespace launcher
{

    class Launcher
    {

        protected:

            /**@brief MPI global communicator*/
            MPI_Comm m_uiCommGlobal;
            
            /**@brief : number of nodes used for a single run*/
            unsigned int m_uiNodesPerJob;

            /**@brief : number of cores per node*/
            unsigned int m_uiCoresPerNode;

            /**@brief : number of jobs to launch*/
            unsigned int m_uiNumJobs;

            MPI_Comm* m_uiCommLocal = NULL;

            unsigned int m_uiLocalSplitIndex=0;


        public:
            
            Launcher(MPI_Comm comm, unsigned int nodes_per_job, unsigned int cores_per_node)
            {
                m_uiCommGlobal   = comm;
                m_uiNodesPerJob  = nodes_per_job;
                m_uiCoresPerNode = cores_per_node;
            }
            
            ~Launcher(){};

            int alloc_sub_communicator(unsigned int njobs)
            {
                int rank_g, npes_g;
                MPI_Comm_rank(m_uiCommGlobal,&rank_g);
                MPI_Comm_size(m_uiCommGlobal,&npes_g);

                m_uiNumJobs=njobs;

                if(npes_g != m_uiNumJobs*m_uiNodesPerJob*m_uiCoresPerNode)
                {
                    if(!rank_g)
                        std::cout<<"GLOBAL CORES : "<<npes_g<<" not enough to run "<<m_uiNumJobs<<" jobs with "<<m_uiNodesPerJob<<" nodes with "<<m_uiCoresPerNode<<" cores per node"<<std::endl;
                    
                    MPI_Abort(m_uiCommGlobal,0);
                }

                const unsigned int cores_per_job = m_uiCoresPerNode * m_uiNodesPerJob;
                unsigned int split_index         = rank_g / cores_per_job;
                m_uiLocalSplitIndex              =  split_index;
                deallocate_sub_communicator();
                m_uiCommLocal = new MPI_Comm();
                MPI_Comm_split(m_uiCommGlobal,split_index,rank_g,m_uiCommLocal);
                return 0;

            }

            int deallocate_sub_communicator()
            {
                if(m_uiCommLocal !=NULL)
                {
                    MPI_Comm_free(m_uiCommLocal);
                    m_uiCommLocal=NULL;
                }

                return 0;

            }

            /**@brief get the subcommunicators list (use with care :) )*/
            inline const MPI_Comm* get_sub_communicator() const { return m_uiCommLocal; }

            /**@brief get the global communicator*/
            inline MPI_Comm get_global_communicator() const { return m_uiCommGlobal; }

            inline unsigned int get_comm_split_index() const { return m_uiLocalSplitIndex;}


    };



}// end of namespcae launcher. 