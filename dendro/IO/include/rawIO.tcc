
    template <typename T>
    void io::varToRawData(const ot::Mesh* pMesh, const T** vars, unsigned int numVars,const char** varNames,const char* fPrefix)
    {

        char fName[256];
        sprintf(fName,"%s.bin",fPrefix);
        std::ofstream ofile(fName,std::ios::binary);

        unsigned int globalNpes = pMesh->getMPICommSizeGlobal();


        const unsigned int eleLocalBegin=pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd=pMesh->getElementLocalEnd();

        const unsigned int Nx=(1u<<m_uiMaxDepth)+1;
        const unsigned int Ny=(1u<<m_uiMaxDepth)+1;
        const unsigned int Nz=(1u<<m_uiMaxDepth)+1;

        const DendroIntL N=Nx*Ny*Nz;

        const unsigned int nPe=pMesh->getNumNodesPerElement();

        const unsigned int nx=pMesh->getElementOrder()+1;
        const unsigned int ny=pMesh->getElementOrder()+1;
        const unsigned int nz=pMesh->getElementOrder()+1;

        const ot::TreeNode* allElements=&(*(pMesh->getAllElements().begin()));

        T* nodalValues=new T[nPe];
        T* regGridBuffer=new T[N];

        ofile.write((char*) &Nx, sizeof(unsigned int));
        ofile.write((char*) &Ny, sizeof(unsigned int));
        ofile.write((char*) &Nz, sizeof(unsigned int));
        ofile.write((char*) &numVars, sizeof(unsigned int));

        unsigned int level,rxI,ryI,rzI;
        unsigned int skipIndex;

        for(unsigned int v=0;v<numVars;v++)
        {
            for(unsigned int ele=eleLocalBegin;ele<eleLocalEnd;ele++)
            {
                level=allElements[ele].getLevel();
                rxI=allElements[ele].minX();
                ryI=allElements[ele].minY();
                rzI=allElements[ele].minZ();

                // get the elemental nodal values.
                pMesh->getElementNodalValues(vars[v],nodalValues,ele);
                skipIndex=(1u<<(m_uiMaxDepth-(level+MAXDEAPTH_LEVEL_DIFF+1)))-1;

                for(unsigned int k=0;k<nz;k++)
                    for(unsigned int j=0;j<ny;j++)
                        for(unsigned int i=0;i<nx;i++)
                        {
                            regGridBuffer[(rzI+k) *(Ny*Nx) + (ryI+j)*(Nx) + (rxI+i) ] = nodalValues[ k*(ny*nx) + j*(ny) + i];
                           /* for(unsigned int sk=0;sk<skipIndex;sk++)
                                for(unsigned int sy=0;sy<skipIndex;sy++)
                                    for(unsigned int si=0;si<skipIndex;si++)
                                    {
                                       regGridBuffer[(rzI+sk) *(Ny*Nx) + (ryI+sy)*(Nx) + (rxI+si) ]= nodalValues[ k*(ny*nx) + j*(ny) + i];
                                    }*/


                        }

            }

            ofile.write((char*)&regGridBuffer[0], sizeof(double)*N);

        }



        delete [] nodalValues;
        delete [] regGridBuffer;

        ofile.close();



    }


    template<typename T>
    void io::dump_array(const  T* const in, unsigned int sz, const char* fPrefix, MPI_Comm comm)
    {

        int rank, npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &npes);

        std::ofstream mfile;
        char fName[256];
        sprintf(fName,"%s_%d_%d.dump",fPrefix,rank,npes);
        mfile.open(fName,std::ios::app);

        if(!mfile)
        {
            std::cout<<"line: "<<__LINE__<<" fopen failed"<<std::endl;
            return ;

        }

        mfile<<"[i]\t val"<<std::endl;
        for(unsigned int i=0; i< sz; i++)
            mfile<<"["<<i<<"]\t val : "<<in[i]<<std::endl;

        mfile.close();
        return;


    }