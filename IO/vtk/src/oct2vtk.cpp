//
// Created by milinda on 5/30/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "oct2vtk.h"

namespace io
{
    namespace vtk
    {
        void mesh2vtk(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,double * filedData,unsigned int numPointdata, const char **pointDataNames, double **pointData)
        {

            if(!(pMesh->isActive())) return ;

            unsigned int rank = pMesh->getMPIRank();
            unsigned int npes = pMesh->getMPICommSize();

            char fname[FNAME_LENGTH];
            sprintf(fname,"%s_%d_%d.vtk",fPrefix,rank,npes);

            fp = fopen(fname,"w+");
            if(fp==NULL) {
                std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                return ;
            }

           const std::vector<ot::TreeNode>& pElements = pMesh->getAllElements();

           write_vtu_header();
           fprintf(fp," DATASET UNSTRUCTURED_GRID\n");

           unsigned int dim            = m_uiDim;
           unsigned int num_vertices   = pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement();
           unsigned int num_cells      = pMesh->getNumLocalMeshElements();
           int      num_cells_elements = num_cells * NUM_CHILDREN + num_cells;
           unsigned int nPe            = pMesh->getNumNodesPerElement();
           unsigned int eleOrder       = pMesh->getElementOrder();
           double sz;

           fprintf(fp,"POINTS %d float\n",num_vertices);

           for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++)
           {
               sz = 1u<<(m_uiMaxDepth-pElements[ele].getLevel());

               for(unsigned int k=0;k<(eleOrder+1);k++)
                 for(unsigned int j=0;j<(eleOrder+1);j++)
                   for(unsigned int i=0;i<(eleOrder+1);i++)
                   {
                       fprintf(fp,"%f %f %f \n",(float)(pElements[ele].getX()+i*(sz/eleOrder)),(float)(pElements[ele].getY()+j*(sz/eleOrder)),(float)(pElements[ele].getZ()+k*(sz/eleOrder)));

                   }


           }

            fprintf(fp,"       \n");
            fprintf(fp, "CELLS %d %d\n", num_cells, num_cells_elements);

            for (unsigned int i = 0 ; i < num_cells ; i++)
            {


                fprintf(fp,"%lld %lld %lld %lld %lld %lld %lld %lld %lld\n",(DendroIntL)NUM_CHILDREN,
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0));



            }

            fprintf(fp,"       \n");
            fprintf(fp, "CELL_TYPES %d\n", num_cells);

            for (int i = 0; i < num_cells; i++) {
                fprintf(fp,"%d\n",VTK_HEXAHEDRON);
            }

            fprintf(fp,"  n");
            fclose(fp);
            fp = NULL;

        }



        void mesh2vtu(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames,const double **pointData,bool isDGPData)
        {


            if(!(pMesh->isActive())) return ;

            unsigned int rank = pMesh->getMPIRank();
            unsigned int npes = pMesh->getMPICommSize();

            unsigned int retval;

            Point dmin = pMesh->getDomainMinPt();
            Point dmax = pMesh->getDomainMaxPt();
            const double invRg = 1.0/( (double) (1u<<m_uiMaxDepth));

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp = fopen(fname,"w+");
            if(fp==NULL) {
                std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                return ;
            }

            const std::vector<ot::TreeNode>& pElements = pMesh->getAllElements();
            const std::vector<unsigned int>& e2N       = pMesh->getE2NMapping();

            unsigned   int dim       = m_uiDim;
            DendroIntL num_vertices  = pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement();
            unsigned   int num_cells = pMesh->getNumLocalMeshElements();

            //int num_cells_elements = num_cells * NUM_CHILDREN + num_cells;
            unsigned int nPe      = pMesh->getNumNodesPerElement();
            unsigned int eleOrder = pMesh->getElementOrder();
            double sz;

            std:: vector<double> nodalVal;
            nodalVal.resize(nPe);

            fprintf(fp,"<?xml version=\"1.0\"?>\n");
            fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
            #ifdef DENDRO_VTU_ZLIB
            fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
            #endif
            fprintf(fp," byte_order=\"LittleEndian\">\n");
            fprintf(fp,"  <UnstructuredGrid>\n");


            fprintf(fp,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);


            if(numFieldData>0 && filedData!=NULL)
            {
                fprintf(fp,"      <FieldData>\n");

                for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                {
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);
                    fprintf(fp,"%f\n",filedData[fdata]);
                    fprintf(fp,"      </DataArray>\n");
                }

                fprintf(fp,"      </FieldData>\n");
            }



            fprintf(fp,"      <Points>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                        #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
            #endif


            #ifdef DENDRO_VTU_ASCII
            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {


                            fprintf(fp, "          %d %d %d\n", VTU_OCT_X_GRID_X(pElements[ele].getX() + i * (sz / eleOrder)),
                                    VTU_OCT_Y_GRID_Y(pElements[ele].getY() + j * (sz / eleOrder)),
                                    VTU_OCT_Z_GRID_Z(pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
            #else

            DENDRO_NODE_COORD_DTYPE* coord_data = new DENDRO_NODE_COORD_DTYPE[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 0] = VTU_OCT_X_GRID_X(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 1] = VTU_OCT_Y_GRID_Y(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 2] = VTU_OCT_Z_GRID_Z(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                  std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std:: endl;
                  fclose(fp);

            }
            delete [] coord_data;
               #endif

            fprintf(fp,"\n");

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Points>\n");
            fprintf(fp,"      <Cells>\n");


               #ifdef  DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
               #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
               #endif




            #ifdef  DENDRO_VTU_ASCII
            for (unsigned int i = 0 ; i < num_cells ; i++)
            {
                fprintf(fp,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                        (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0));
            }

            #else
            DendroIntL * locidx_data = new DendroIntL [num_cells*NUM_CHILDREN];
            for (unsigned int i = 0 ; i < num_cells ; i++)
            {
                locidx_data[i*NUM_CHILDREN+0] = (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+1] = (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+2] = (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+3] = (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+4] = (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+5] = (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+6] = (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+7] = (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);

            }
            retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
                fclose(fp);

            }
            delete [] locidx_data;

            #endif



            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                fprintf(fp," %lld",(NUM_CHILDREN*(il)));
                if (!(sk % 8) && il != num_cells)
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
            DendroIntL * loc_offset = new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1] = NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_offset;


            #endif

            fprintf(fp,"\n        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                fprintf (fp, " %d",VTK_HEXAHEDRON);
                if (!(sk % 20) && il != (num_cells - 1))
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);
            uint8_t * loc_type = new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il] = VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_type;


               #endif
            fprintf (fp, "\n");
            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);

            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il] = rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_rank;


            #endif
            fprintf (fp, "\n");
            fprintf(fp,"        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pElements[il+pMesh->getElementLocalBegin()].getLevel());
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il] = pElements[il+pMesh->getElementLocalBegin()].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_level;

            #endif

            fprintf (fp, "\n");
            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </CellData>\n");


            // writing point data
            if(numPointdata>0 && pointData!=NULL)
            {
                fprintf(fp,"      <PointData>\n");

                #ifndef DENDRO_VTU_ASCII
                double * nodalVal_binary = new double [num_cells*nPe];
                #endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

                #ifdef DENDRO_VTU_ASCII
                  fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                #else
                  fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                #endif


                  const double *tmp_var = pointData[pdata];
                  // note that size of tmp_var should be actual number of nodes in the mesh.
                  //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
                  #ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPData);
                        for(unsigned int node=0;node<nPe;node++)
                            fprintf(fp,"%f ",nodalVal[node]);
                    }
                    fprintf(fp,"\n");
                  #else

                  for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPData);
                        for(unsigned int node=0;node<nPe;node++)
                            nodalVal_binary[(il-pMesh->getElementLocalBegin())*nPe+node] = nodalVal[node];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_cells*nPe);

                    if (retval) {

                        std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std:: endl;
                        fclose(fp);
                    }


                #endif
                    fprintf (fp, "\n");
                    fprintf(fp,"        </DataArray>\n");
                
                }
                fprintf(fp,"      </PointData>\n");


                #ifndef DENDRO_VTU_ASCII
                 delete [] nodalVal_binary;
                #endif
            }




            fprintf(fp,"</Piece>\n");
            fprintf(fp,"</UnstructuredGrid>\n");
            fprintf(fp,"</VTKFile>\n");

            fclose(fp);
            fp = NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp = fopen(fname,"w+");
                if(fp==NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std:: endl;
                    return ;
                }


                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
                #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
                #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");

                fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

                if(numFieldData>0 && filedData!=NULL)
                {
                    fprintf(fp,"    <PFieldData>\n");
                    for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                    {
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);

                    }
                    fprintf(fp,"    </PFieldData>\n");
                }

                fprintf(fp,"    <PPoints>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
                #endif



                fprintf(fp,"    </PPoints>\n");
                fprintf(fp,"    <PCellData>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif



                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif

                fprintf(fp,"    </PCellData>\n");


                // writing point data
                if(numPointdata>0 && pointData!=NULL)
                {
                    fprintf(fp,"      <PPointData>\n");

                    for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                    {

                #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                #else
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                #endif

                    }

                    fprintf(fp,"      </PPointData>\n");

                }
                std::string vtuName = getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp = NULL;




            }


        }



        void mesh2vtuCoarse(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData, bool isDGPData)
        {

            if(!(pMesh->isActive())) return ;

            unsigned int rank = pMesh->getMPIRank();
            unsigned int npes = pMesh->getMPICommSize();

            Point dmin = pMesh->getDomainMinPt();
            Point dmax = pMesh->getDomainMaxPt();
            const double invRg = 1.0/( (double) (1u<<m_uiMaxDepth));

            int retval;

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp = fopen(fname,"w+");
            if(fp==NULL) {
                std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                return ;
            }

            const std::vector<ot::TreeNode>& pElements = pMesh->getAllElements();
            const std::vector<unsigned int>& e2N       = pMesh->getE2NMapping();

            unsigned   int dim       = m_uiDim;
            DendroIntL num_vertices  = pMesh->getNumLocalMeshElements()*NUM_CHILDREN;
            unsigned   int num_cells = pMesh->getNumLocalMeshElements();
            unsigned   int nPe       = pMesh->getNumNodesPerElement();
            unsigned   int eleOrder  = pMesh->getElementOrder();
            double sz;

            std:: vector<double> nodalVal;
            nodalVal.resize(nPe);

            const unsigned int nx = 2;
            const unsigned int ny = 2;
            const unsigned int nz = 2;

            fprintf(fp,"<?xml version=\"1.0\"?>\n");
            fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
            
            #ifdef DENDRO_VTU_ZLIB
            fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
            #endif
            fprintf(fp," byte_order=\"LittleEndian\">\n");
            fprintf(fp,"  <UnstructuredGrid>\n");


            fprintf(fp,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);

            if(numFieldData>0 && filedData!=NULL)
            {
                fprintf(fp,"      <FieldData>\n");

                for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                {
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);
                    fprintf(fp,"%f\n",filedData[fdata]);
                    fprintf(fp,"      </DataArray>\n");
                }

                fprintf(fp,"      </FieldData>\n");
            }



            fprintf(fp,"      <Points>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
            #endif


            fprintf(fp,"          ");

            #ifdef DENDRO_VTU_ASCII
            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k+=eleOrder)
                    for (unsigned int j = 0; j < (eleOrder + 1); j+=eleOrder)
                        for (unsigned int i = 0; i < (eleOrder + 1); i+=eleOrder) {


                            fprintf(fp, "          %d %d %d\n", VTU_OCT_X_GRID_X(pElements[ele].getX() + i * (sz / eleOrder)),
                                    VTU_OCT_Y_GRID_Y(pElements[ele].getY() + j * (sz / eleOrder)),
                                    VTU_OCT_Z_GRID_Z(pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
            #else
            
            DENDRO_NODE_COORD_DTYPE* coord_data = new DENDRO_NODE_COORD_DTYPE[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k+=eleOrder)
                    for (unsigned int j = 0; j < (eleOrder + 1); j+=eleOrder)
                        for (unsigned int i = 0; i < (eleOrder + 1); i+=eleOrder) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 0] = VTU_OCT_X_GRID_X(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 1] = VTU_OCT_Y_GRID_Y(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 2] = VTU_OCT_Z_GRID_Z(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std:: endl;
                fclose(fp);

            }
            delete [] coord_data;
            #endif

            fprintf(fp,"\n");

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Points>\n");
            fprintf(fp,"      <Cells>\n");


            #ifdef  DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
            #endif



            for (unsigned int i = 0 ; i < num_cells ; i++)
            {


                #ifdef  DENDRO_VTU_ASCII
                fprintf(fp,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0));

               #else


                DendroIntL * locidx_data = new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int i = 0 ; i < num_cells ; i++)
                {
                    locidx_data[i*NUM_CHILDREN+0] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+1] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+2] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+3] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0);
                    locidx_data[i*NUM_CHILDREN+4] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+5] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+6] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+7] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0);

                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
                    fclose(fp);

                }
                delete [] locidx_data;


               #endif

            }

            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                fprintf(fp," %lld",(NUM_CHILDREN*il));
                if (!(sk % 8) && il != num_cells)
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);

            DendroIntL * loc_offset = new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1] = NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_offset;

            #endif



            fprintf(fp,"\n        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                fprintf (fp, " %d",VTK_HEXAHEDRON);
                if (!(sk % 20) && il != (num_cells - 1))
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);
            uint8_t * loc_type = new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il] = VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_type;

            #endif

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il] = rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_rank;
            #endif

            fprintf(fp,"        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pElements[il+pMesh->getElementLocalBegin()].getLevel());
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il] = pElements[il+pMesh->getElementLocalBegin()].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_level;
            #endif


            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </CellData>\n");




            // writing point data
            if(numPointdata>0 && pointData!=NULL)
            {
                fprintf(fp,"      <PointData>\n");

                #ifndef DENDRO_VTU_ASCII
                double * nodalVal_binary = new double [num_vertices];
                #endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

                    #ifdef DENDRO_VTU_ASCII
                     fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                    #else
                     fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                    #endif


                    const double *tmp_var = pointData[pdata];
                    // note that size of tmp_var should be actual number of nodes in the mesh.
                    //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
                    #ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var, &(*(nodalVal.begin())), il,isDGPData);
                        for (unsigned int k = 0; k < (eleOrder + 1); k += eleOrder)
                            for (unsigned int j = 0; j < (eleOrder + 1); j += eleOrder)
                                for (unsigned int i = 0; i < (eleOrder + 1); i += eleOrder)
                                    fprintf(fp,"%f ",nodalVal[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]);
                    }
                    fprintf(fp,"\n");
                    #else



                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPData);
                        for (unsigned int k = 0; k < (eleOrder + 1); k += eleOrder)
                            for (unsigned int j = 0; j < (eleOrder + 1); j += eleOrder)
                                for (unsigned int i = 0; i < (eleOrder + 1); i += eleOrder)
                                    nodalVal_binary[(il-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder)] = nodalVal[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_vertices);

                    if (retval) {

                        std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std:: endl;
                        fclose(fp);
                    }


                    #endif
                    fprintf (fp, "\n");
                    fprintf(fp,"        </DataArray>\n");
                }
                fprintf(fp,"      </PointData>\n");


                #ifndef DENDRO_VTU_ASCII
                delete [] nodalVal_binary;
                #endif
            }



            fprintf(fp,"</Piece>\n");
            fprintf(fp,"</UnstructuredGrid>\n");
            fprintf(fp,"</VTKFile>\n");

            fclose(fp);
            fp = NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp = fopen(fname,"w+");
                if(fp==NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std:: endl;
                    return ;
                }


                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
                #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
                #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");

                fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

                if(numFieldData>0 && filedData!=NULL)
                {
                    fprintf(fp,"    <PFieldData>\n");
                    for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                    {
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);

                    }
                    fprintf(fp,"    </PFieldData>\n");
                }

                fprintf(fp,"    <PPoints>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
                #endif



                fprintf(fp,"    </PPoints>\n");
                fprintf(fp,"    <PCellData>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif



                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif

                fprintf(fp,"    </PCellData>\n");


                // writing point data
                if(numPointdata>0 && pointData!=NULL)
                {
                    fprintf(fp,"      <PPointData>\n");

                    for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                    {

                #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                #else
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                #endif

                    }

                    fprintf(fp,"      </PPointData>\n");

                }
                std::string vtuName = getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp = NULL;




            }


        }


        void mesh2vtuFine(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData, unsigned int nCellData, const char** cellDNames, const double** cellData, bool isDGPdata)
        {

            if(!(pMesh->isActive())) return ;

            unsigned int rank = pMesh->getMPIRank();
            unsigned int npes = pMesh->getMPICommSize();

            Point dmin = pMesh->getDomainMinPt();
            Point dmax = pMesh->getDomainMaxPt();
            const double invRg = 1.0/( (double) (1u<<m_uiMaxDepth));

            int retval;

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp = fopen(fname,"w+");
            if(fp==NULL) {
                std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                return ;
            }

            const std::vector<ot::TreeNode>& pElements = pMesh->getAllElements();
            const std::vector<unsigned int>& e2N       = pMesh->getE2NMapping();

            // note this works only for the 3D grids.
            const unsigned int dim        = m_uiDim;
            const unsigned int nPe        = pMesh->getNumNodesPerElement();
            const unsigned int eleOrder   = pMesh->getElementOrder();
            const unsigned int ePe        = eleOrder*eleOrder*eleOrder;
            const unsigned int num_cells  = pMesh->getNumLocalMeshElements()*ePe;
            const DendroIntL num_vertices = pMesh->getNumLocalMeshElements()*nPe;
            double sz;

            std:: vector<double> nodalVal;
            nodalVal.resize(nPe);

            fprintf(fp,"<?xml version=\"1.0\"?>\n");
            fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
            
            #ifdef DENDRO_VTU_ZLIB
            fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
            #endif
            fprintf(fp," byte_order=\"LittleEndian\">\n");
            fprintf(fp,"  <UnstructuredGrid>\n");


            fprintf(fp,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);


            if(numFieldData>0 && filedData!=NULL)
            {
                fprintf(fp,"      <FieldData>\n");

                for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                {
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);
                    fprintf(fp,"%f\n",filedData[fdata]);
                    fprintf(fp,"      </DataArray>\n");
                }

                fprintf(fp,"      </FieldData>\n");
            }



            fprintf(fp,"      <Points>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
            #endif


            #ifdef DENDRO_VTU_ASCII
            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {


                            fprintf(fp, "          %d %d %d\n", VTU_OCT_X_GRID_X(pElements[ele].getX() + i * (sz / eleOrder)),
                                    VTU_OCT_Y_GRID_Y(pElements[ele].getY() + j * (sz / eleOrder)),
                                    VTU_OCT_Z_GRID_Z(pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
            #else

            DENDRO_NODE_COORD_DTYPE* coord_data = new DENDRO_NODE_COORD_DTYPE[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 0] = VTU_OCT_X_GRID_X(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 1] = VTU_OCT_Y_GRID_Y(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 2] = VTU_OCT_Z_GRID_Z(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std:: endl;
                fclose(fp);

            }
            delete [] coord_data;
            #endif


            fprintf(fp,"\n");

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Points>\n");
            fprintf(fp,"      <Cells>\n");



            #ifdef  DENDRO_VTU_ASCII

            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);

            for (unsigned int ele = 0 ; ele < pMesh->getNumLocalMeshElements() ; ele++) {
                for(unsigned int ek=0;ek<(eleOrder);ek++)
                    for(unsigned int ej=0;ej<(eleOrder);ej++)
                        for(unsigned int ei=0;ei<(eleOrder);ei++)
                        {
                            fprintf(fp,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                                    (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei),
                                    (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1)),
                                    (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1)),
                                    (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei),
                                    (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei),
                                    (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1)),
                                    (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1)),
                                    (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei));

                        }


             }

            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);

                DendroIntL * locidx_data = new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int ele = 0 ; ele < pMesh->getNumLocalMeshElements() ; ele++)
                {
                    for(unsigned int ek=0;ek<(eleOrder);ek++)
                        for(unsigned int ej=0;ej<(eleOrder);ej++)
                            for(unsigned int ei=0;ei<(eleOrder);ei++)
                            {

                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+0] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+1] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+2] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+3] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+4] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+5] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+6] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+7] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                            }




                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
                    fclose(fp);

                }
                delete [] locidx_data;

            #endif




            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                fprintf(fp," %lld",(NUM_CHILDREN*il));
                if (!(sk % 8) && il != num_cells)
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
            DendroIntL * loc_offset = new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1] = NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_offset;
            #endif



            fprintf(fp,"\n        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                fprintf (fp, " %d",VTK_HEXAHEDRON);
                if (!(sk % 20) && il != (num_cells - 1))
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);
            uint8_t * loc_type = new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il] = VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_type;
            
            #endif
            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

            #ifndef DENDRO_VTU_ASCII
                unsigned int * cell_tmp = new unsigned int [num_cells];
                double * cell_dtmp =NULL;
                
                if(nCellData>0 && cellData!=NULL) 
                    cell_dtmp=new double[num_cells];
            #endif

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            
            for (unsigned int il = 0; il < num_cells; ++il)
                cell_tmp[il] = rank;

            retval = vtk_write_binary (fp, (char *) cell_tmp, sizeof (*cell_tmp)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std:: endl;
                fclose(fp);
            }
            #endif
            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
               for(unsigned int w=0;w<ePe;w++)
                  fprintf(fp,"%d ",(pElements[il].getLevel()));
            }

            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
            for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                for(unsigned int w=0;w<ePe;w++)
                    cell_tmp[(il-pMesh->getElementLocalBegin())*ePe+w] = pElements[il].getLevel();
            }

            retval = vtk_write_binary (fp, (char *) cell_tmp, sizeof (*cell_tmp)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                fclose(fp);
            }
            #endif

            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");

            if(nCellData > 0 && cellData!=NULL) 
            {

                for(unsigned int v=0; v < nCellData; v++)
                {

                    #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE, cellDNames[v], DENDRO_FORMAT_ASCII);
                        fprintf(fp,"         ");
                        for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        for(unsigned int w=0;w<ePe;w++)
                            fprintf(fp,"%f ", cellData[v][il] );
                        }

                        fprintf(fp,"\n");
                    #else
                        fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,cellDNames[v],DENDRO_FORMAT_BINARY);

                        for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                            for(unsigned int w=0;w<ePe;w++)    
                             cell_dtmp[(il-pMesh->getElementLocalBegin())*ePe+w] = cellData[v][il];
                        }

                        retval = vtk_write_binary (fp, (char *) (cell_dtmp), sizeof (*cell_dtmp)*num_cells);
                        if (retval) {
                            std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                            fclose(fp);
                        }

                    #endif

                    fprintf(fp,"\n");
                    fprintf(fp,"        </DataArray>\n");

                }

            }


            #ifndef DENDRO_VTU_ASCII
                delete [] cell_tmp;
                
                if(nCellData>0 && cellData!=NULL) 
                    delete [] cell_dtmp;
            #endif

            
            fprintf(fp,"      </CellData>\n");




            if(numPointdata>0 && pointData!=NULL)
            {
                fprintf(fp,"      <PointData>\n");

                #ifndef DENDRO_VTU_ASCII
                double * nodalVal_binary = new double [num_vertices];
                #endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

                    #ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                    #else
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                    #endif


                    const double *tmp_var = pointData[pdata];
                    // note that size of tmp_var should be actual number of nodes in the mesh.
                    //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
                    #ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPdata);
                        for(unsigned int node=0;node<nPe;node++)
                            fprintf(fp,"%f ",nodalVal[node]);
                    }
                    fprintf(fp,"\n");
                    #else



                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPdata);
                        for(unsigned int node=0;node<nPe;node++)
                            nodalVal_binary[(il-pMesh->getElementLocalBegin())*nPe+node] = nodalVal[node];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_vertices);

                    if (retval) {

                        std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std:: endl;
                        fclose(fp);
                    }


                    #endif
                    fprintf (fp, "\n");
                    fprintf(fp,"        </DataArray>\n");
                }
                fprintf(fp,"      </PointData>\n");


                #ifndef DENDRO_VTU_ASCII
                delete [] nodalVal_binary;
                #endif
            }

            fprintf (fp, "\n");
            fprintf(fp,"</Piece>\n");
            fprintf(fp,"</UnstructuredGrid>\n");
            fprintf(fp,"</VTKFile>\n");

            fclose(fp);
            fp = NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp = fopen(fname,"w+");
                if(fp==NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std:: endl;
                    return ;
                }


                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
                #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
                #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");

                fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

                if(numFieldData>0 && filedData!=NULL)
                {
                    fprintf(fp,"    <PFieldData>\n");
                    for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                    {
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);

                    }
                    fprintf(fp,"    </PFieldData>\n");
                }

                fprintf(fp,"    <PPoints>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
                #endif



                fprintf(fp,"    </PPoints>\n");
                fprintf(fp,"    <PCellData>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif



                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif

                if(nCellData>0 && cellData!=NULL)
                {
                    
                    for(unsigned int v=0; v < nCellData; v++)
                    {
                        #ifdef DENDRO_VTU_ASCII
                            fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,cellDNames[v],DENDRO_FORMAT_ASCII);
                        #else
                            fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,cellDNames[v],DENDRO_FORMAT_BINARY);
                        #endif
                    }

                }

                fprintf(fp,"    </PCellData>\n");


                // writing point data
                if(numPointdata>0 && pointData!=NULL)
                {
                    fprintf(fp,"      <PPointData>\n");

                    for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                    {

                        #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                        #else
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                        #endif

                    }

                    fprintf(fp,"      </PPointData>\n");

                }

                


                std::string vtuName = getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp = NULL;




            }
        }


        void oct2vtu(const ot::TreeNode *pNodes,const unsigned int nSize,const char* fPrefix, MPI_Comm comm)
        {
            int rank,npes;
            MPI_Comm_rank(comm,&rank);
            MPI_Comm_size(comm,&npes);

            unsigned int retval;
            unsigned int eleOrder = 1;


            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp = fopen(fname,"w+");
            if(fp==NULL) {
                std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                return ;
            }


            unsigned   int dim       = m_uiDim;
            DendroIntL num_vertices  = nSize*NUM_CHILDREN;
            unsigned   int num_cells = nSize;
            double sz;

            fprintf(fp,"<?xml version=\"1.0\"?>\n");
            fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
            
            #ifdef DENDRO_VTU_ZLIB
            fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
            #endif
            
            fprintf(fp," byte_order=\"LittleEndian\">\n");
            fprintf(fp,"  <UnstructuredGrid>\n");


            fprintf(fp,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);

            fprintf(fp,"      <Points>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
            #endif


            fprintf(fp,"          ");

            const unsigned int nx = 2;
            const unsigned int ny = 2;
            const unsigned int nz = 2;


            #ifdef DENDRO_VTU_ASCII
            for(unsigned int ele=0;ele<nSize;ele++) {
                sz = 1u << (m_uiMaxDepth - pNodes[ele].getLevel());
                for(unsigned int k=0;k<nz;k++)
                   for(unsigned int j=0;j<ny;j++)
                     for(unsigned int i=0;i<nx;i++) {
                            fprintf(fp, "          %d %d %d\n", (pNodes[ele].getX() + i * sz), (pNodes[ele].getY() + j * sz), (pNodes[ele].getZ() + k * sz));
                        }
            }
            #else

            unsigned int* coord_data = new unsigned int[num_vertices*m_uiDim];


            for(unsigned int ele=0;ele<nSize;ele++) {

                sz = 1u << (m_uiMaxDepth - pNodes[ele].getLevel());

                for(unsigned int k=0;k<nz;k++)
                   for(unsigned int j=0;j<ny;j++)
                     for(unsigned int i=0;i<nx;i++)
                     {
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 0] = (pNodes[ele].getX()+i*sz);
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 1] = (pNodes[ele].getY()+j*sz);
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 2] = (pNodes[ele].getZ()+k*sz);
                     }


            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std:: endl;
                fclose(fp);

            }
            delete [] coord_data;
            #endif

            fprintf(fp,"\n");

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Points>\n");
            fprintf(fp,"      <Cells>\n");


            #ifdef  DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
            #endif



            for (unsigned int i = 0 ; i < num_cells ; i++)
            {


               #ifdef  DENDRO_VTU_ASCII
                fprintf(fp,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1),
                        (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0));

               #else


                DendroIntL * locidx_data = new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int i = 0 ; i < num_cells ; i++)
                {
                    locidx_data[i*NUM_CHILDREN+0] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+1] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+2] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+3] = (DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0);
                    locidx_data[i*NUM_CHILDREN+4] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+5] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+6] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+7] = (DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0);

                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
                    fclose(fp);

                }
                delete [] locidx_data;


               #endif

            }

            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                fprintf(fp," %lld",(NUM_CHILDREN*il));
                if (!(sk % 8) && il != num_cells)
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);

            DendroIntL * loc_offset = new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1] = NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_offset;

            #endif



            fprintf(fp,"\n        </DataArray>\n");



            #ifdef DENDRO_VTU_ASCII
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                fprintf (fp, " %d",VTK_HEXAHEDRON);
                if (!(sk % 20) && il != (num_cells - 1))
                    fprintf(fp,"\n         ");
            }
            fprintf(fp,"\n");
            #else
            fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);
            uint8_t * loc_type = new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il] = VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_type;

            #endif

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il] = rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_rank;
            #endif

            fprintf(fp,"        </DataArray>\n");


            #ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pNodes[il].getLevel());
            }
            fprintf(fp,"\n");
            #else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level = new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il] = pNodes[il].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                fclose(fp);
            }
            delete [] loc_level;
            #endif


            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </CellData>\n");
            fprintf(fp,"</Piece>\n");
            fprintf(fp,"</UnstructuredGrid>\n");
            fprintf(fp,"</VTKFile>\n");

            fclose(fp);
            fp = NULL;

            if(!rank) {
                sprintf(fname, "%s.pvtu", fPrefix);

                fp = fopen(fname, "w+");
                if (fp == NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std:: endl;
                    return;
                }


                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
               #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
               #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");

                fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");


                fprintf(fp,"    <PPoints>\n");

               #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
               #else
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
               #endif



                fprintf(fp,"    </PPoints>\n");
                fprintf(fp,"    <PCellData>\n");

               #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
               #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
               #endif



               #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
               #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
               #endif

                fprintf(fp,"    </PCellData>\n");
                std::string vtuName = getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp = NULL;


            }



        }


        void mesh2vtu_slice(const ot::Mesh *pMesh, unsigned int s_val[], unsigned int s_normal[], const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData,bool isDGPData)
        {
             if(!(pMesh->isActive())) return ;

            unsigned int rank = pMesh->getMPIRank();
            unsigned int npes = pMesh->getMPICommSize();
            MPI_Comm activeComm = pMesh->getMPICommunicator();

            Point dmin = pMesh->getDomainMinPt();
            Point dmax = pMesh->getDomainMaxPt();
            const double invRg = 1.0/( (double) (1u<<m_uiMaxDepth));

            int retval;

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            std::vector<unsigned int> sids;
            ot::slice_mesh(pMesh,s_val,s_normal,sids);

            bool sid_empty = sids.empty();
            bool* sid_empty_g = new bool[npes];

            
           
            
            if(!sid_empty)
            {

                fp = fopen(fname,"w+");
                if(fp==NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std:: endl;
                    return ;
                }

                const std::vector<ot::TreeNode>& pElements = pMesh->getAllElements();
                const std::vector<unsigned int>& e2N       = pMesh->getE2NMapping();

                // note this works only for the 3D grids.
                const unsigned int dim        = m_uiDim;
                const unsigned int nPe        = pMesh->getNumNodesPerElement();
                const unsigned int eleOrder   = pMesh->getElementOrder();
                const unsigned int ePe        = eleOrder*eleOrder*eleOrder;
                const unsigned int num_cells  = sids.size()*ePe;
                const DendroIntL num_vertices = sids.size()*nPe;
                double sz;

                std:: vector<double> nodalVal;
                nodalVal.resize(nPe);

            

                

                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
                
                #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
                #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");
                fprintf(fp,"  <UnstructuredGrid>\n");


                fprintf(fp,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);


                if(numFieldData>0 && filedData!=NULL)
                {
                    fprintf(fp,"      <FieldData>\n");

                    for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                    {
                        fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\" format=\"%s\">\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);
                        fprintf(fp,"%f\n",filedData[fdata]);
                        fprintf(fp,"      </DataArray>\n");
                    }

                    fprintf(fp,"      </FieldData>\n");
                }



                fprintf(fp,"      <Points>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
                #endif


                #ifdef DENDRO_VTU_ASCII
                for(unsigned int s = 0 ; s < sids.size(); s++) {
                    const unsigned int ele = sids[s];
                    sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                    for (unsigned int k = 0; k < (eleOrder + 1); k++)
                        for (unsigned int j = 0; j < (eleOrder + 1); j++)
                            for (unsigned int i = 0; i < (eleOrder + 1); i++) {


                                fprintf(fp, "          %d %d %d\n", VTU_OCT_X_GRID_X(pElements[ele].getX() + i * (sz / eleOrder)),
                                        VTU_OCT_Y_GRID_Y(pElements[ele].getY() + j * (sz / eleOrder)),
                                        VTU_OCT_Z_GRID_Z(pElements[ele].getZ() + k * (sz / eleOrder)));

                            }
                }
                #else

                DENDRO_NODE_COORD_DTYPE* coord_data = new DENDRO_NODE_COORD_DTYPE[num_vertices*m_uiDim];

                for(unsigned int s = 0 ; s < sids.size(); s++) {
                    const unsigned int ele = sids[s];
                    sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                    for (unsigned int k = 0; k < (eleOrder + 1); k++)
                        for (unsigned int j = 0; j < (eleOrder + 1); j++)
                            for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                                coord_data[((s)*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 0] = VTU_OCT_X_GRID_X(pElements[ele].getX()+i*(sz/eleOrder));
                                coord_data[((s)*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 1] = VTU_OCT_Y_GRID_Y(pElements[ele].getY()+j*(sz/eleOrder));
                                coord_data[((s)*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 2] = VTU_OCT_Z_GRID_Z(pElements[ele].getZ()+k*(sz/eleOrder));
                            }
                }

                retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std:: endl;
                    fclose(fp);

                }
                delete [] coord_data;
                #endif


                fprintf(fp,"\n");

                fprintf(fp,"        </DataArray>\n");
                fprintf(fp,"      </Points>\n");
                fprintf(fp,"      <Cells>\n");



                #ifdef  DENDRO_VTU_ASCII

                fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);

                for (unsigned int ele = 0 ; ele < sids.size() ; ele++) {
                    for(unsigned int ek=0;ek<(eleOrder);ek++)
                        for(unsigned int ej=0;ej<(eleOrder);ej++)
                            for(unsigned int ei=0;ei<(eleOrder);ei++)
                            {
                                fprintf(fp,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                                        (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei),
                                        (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1)),
                                        (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1)),
                                        (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei),
                                        (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei),
                                        (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1)),
                                        (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1)),
                                        (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei));

                            }


                }

                #else
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);

                    DendroIntL * locidx_data = new DendroIntL [num_cells*NUM_CHILDREN];
                    for (unsigned int ele = 0 ; ele < sids.size() ; ele++)
                    {
                        for(unsigned int ek=0;ek<(eleOrder);ek++)
                            for(unsigned int ej=0;ej<(eleOrder);ej++)
                                for(unsigned int ei=0;ei<(eleOrder);ei++)
                                {
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+0] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+1] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+2] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+3] = (DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+4] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+5] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+6] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                    locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+7] = (DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                                }

                    }

                    retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                    if (retval) {

                        std:: cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std:: endl;
                        fclose(fp);

                    }
                    delete [] locidx_data;

                #endif




                fprintf(fp,"\n");
                fprintf(fp,"        </DataArray>\n");



                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
                fprintf(fp,"         ");
                for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                    fprintf(fp," %lld",(NUM_CHILDREN*il));
                    if (!(sk % 8) && il != num_cells)
                        fprintf(fp,"\n         ");
                }
                fprintf(fp,"\n");
                
                #else
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
                DendroIntL * loc_offset = new DendroIntL [num_cells];
                for (unsigned int il = 1; il <= (num_cells); il++)
                    loc_offset[il-1] = NUM_CHILDREN*il;

                retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std:: endl;
                    fclose(fp);
                }
                delete [] loc_offset;
                #endif



                fprintf(fp,"\n        </DataArray>\n");


                #ifdef DENDRO_VTU_ASCII
                fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
                fprintf(fp,"         ");
                for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                    fprintf (fp, " %d",VTK_HEXAHEDRON);
                    if (!(sk % 20) && il != (num_cells - 1))
                        fprintf(fp,"\n         ");
                }
                fprintf(fp,"\n");
                #else
                fprintf (fp, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);
                uint8_t * loc_type = new uint8_t [num_cells];
                for (unsigned int il = 0; il < num_cells; ++il)
                    loc_type[il] = VTK_HEXAHEDRON;

                retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std:: endl;
                    fclose(fp);
                }
                delete [] loc_type;
                
                #endif
                fprintf(fp,"\n");
                fprintf(fp,"        </DataArray>\n");
                fprintf(fp,"      </Cells>\n");

                // writing cell data
                fprintf(fp,"      <CellData>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                fprintf (fp,"         ");
                for (unsigned int il = 0; il < num_cells; ++il) {
                    fprintf(fp,"%d ",rank);
                }
                fprintf(fp,"\n");
                #else
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);

                unsigned int * loc_rank = new unsigned int [num_cells];
                for (unsigned int il = 0; il < num_cells; ++il)
                    loc_rank[il] = rank;

                retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std:: endl;
                    fclose(fp);
                }
                delete [] loc_rank;
                #endif
                fprintf(fp,"\n");
                fprintf(fp,"        </DataArray>\n");


                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                fprintf(fp,"         ");
                for(unsigned int s = 0 ; s < sids.size(); s++) {
                const unsigned int il = sids[s];
                for(unsigned int w=0;w<ePe;w++)
                    fprintf(fp,"%d ",(pElements[il].getLevel() + (eleOrder>>1u)));
                }

                fprintf(fp,"\n");
                #else
                fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                unsigned int * loc_level = new unsigned int [num_cells];
                for(unsigned int s = 0 ; s < sids.size(); s++) {
                    const unsigned int il = sids[s];
                    for(unsigned int w=0;w<ePe;w++)
                        loc_level[(s)*ePe+w] = pElements[il].getLevel();
                }

                retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
                if (retval) {

                    std:: cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std:: endl;
                    fclose(fp);
                }
                delete [] loc_level;
                #endif

                fprintf(fp,"\n");
                fprintf(fp,"        </DataArray>\n");
                fprintf(fp,"      </CellData>\n");




                if(numPointdata>0 && pointData!=NULL)
                {
                    fprintf(fp,"      <PointData>\n");

                    #ifndef DENDRO_VTU_ASCII
                    double * nodalVal_binary = new double [num_vertices];
                    #endif

                    for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                    {

                        #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                        #else
                        fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                        #endif


                        const double *tmp_var = pointData[pdata];
                        // note that size of tmp_var should be actual number of nodes in the mesh.
                        //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
                        #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"         ");
                        for(unsigned int s = 0 ; s < sids.size(); s++) {
                            const unsigned int il = sids[s];
                            pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPData);
                            for(unsigned int node=0;node<nPe;node++)
                                fprintf(fp,"%f ",nodalVal[node]);
                        }
                        fprintf(fp,"\n");
                        #else



                        for(unsigned int s = 0 ; s < sids.size(); s++) {
                            const unsigned int il = sids[s];
                            pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il,isDGPData);
                            for(unsigned int node=0;node<nPe;node++)
                                nodalVal_binary[(s)*nPe+node] = nodalVal[node];
                        }


                        retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_vertices);

                        if (retval) {

                            std:: cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std:: endl;
                            fclose(fp);
                        }


                        #endif
                        fprintf (fp, "\n");
                        fprintf(fp,"        </DataArray>\n");
                    }
                    fprintf(fp,"      </PointData>\n");


                    #ifndef DENDRO_VTU_ASCII
                    delete [] nodalVal_binary;
                    #endif
                }

                fprintf (fp, "\n");
                fprintf(fp,"</Piece>\n");
                fprintf(fp,"</UnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");

                fclose(fp);
                fp = NULL;

            }
            
            
            par::Mpi_Gather(&sid_empty,sid_empty_g,1,0,activeComm);
            
            

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp = fopen(fname,"w+");
                if(fp==NULL) {
                    std:: cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std:: endl;
                    return ;
                }


                fprintf(fp,"<?xml version=\"1.0\"?>\n");
                fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
                #ifdef DENDRO_VTU_ZLIB
                fprintf(fp," compressor=\"vtkZLibDataCompressor\"");
                #endif
                fprintf(fp," byte_order=\"LittleEndian\">\n");

                fprintf(fp,"  <PUnstructuredGrid GhostLevel=\"0\">\n");

                if(numFieldData>0 && filedData!=NULL)
                {
                    fprintf(fp,"    <PFieldData>\n");
                    for(unsigned int fdata=0;fdata<numFieldData;fdata++)
                    {
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",DENDRO_NODE_VAR_FLOAT,filedDataNames[fdata],DENDRO_FORMAT_ASCII);

                    }
                    fprintf(fp,"    </PFieldData>\n");
                }

                fprintf(fp,"    <PPoints>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
                #endif



                fprintf(fp,"    </PPoints>\n");
                fprintf(fp,"    <PCellData>\n");

                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif



                #ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_ASCII);
                #else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_INT,DENDRO_FORMAT_BINARY);
                #endif

                fprintf(fp,"    </PCellData>\n");


                // writing point data
                if(numPointdata>0 && pointData!=NULL)
                {
                    fprintf(fp,"      <PPointData>\n");

                    for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                    {

                        #ifdef DENDRO_VTU_ASCII
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
                        #else
                        fprintf(fp,"        <PDataArray type=\"%s\" Name=\"%s\"" " format=\"%s\"/>\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
                        #endif

                    }

                    fprintf(fp,"      </PPointData>\n");

                }
                std::string vtuName = getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    if(!sid_empty_g[proc])
                        fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp = NULL;




            }

            delete [] sid_empty_g;





        }
        
        

    } // end of namespace vtk

} // end of namespace io
