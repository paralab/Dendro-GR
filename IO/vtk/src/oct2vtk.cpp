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

            unsigned int rank=pMesh->getMPIRank();
            unsigned int npes=pMesh->getMPICommSize();

            char fname[FNAME_LENGTH];
            sprintf(fname,"%s_%d_%d.vtk",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }

           const std::vector<ot::TreeNode> pElements=pMesh->getAllElements();

           write_vtu_header();
           fprintf(fp," DATASET UNSTRUCTURED_GRID\n");

           unsigned int dim=m_uiDim;
           unsigned int num_vertices=pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement();
           unsigned int num_cells=pMesh->getNumLocalMeshElements();
           int num_cells_elements = num_cells * NUM_CHILDREN + num_cells;
           unsigned int nPe=pMesh->getNumNodesPerElement();
           unsigned int eleOrder=pMesh->getElementOrder();
           unsigned int sz;

           fprintf(fp,"POINTS %d float\n",num_vertices);

           for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++)
           {
               sz=1u<<(m_uiMaxDepth-pElements[ele].getLevel());

               for(unsigned int k=0;k<(eleOrder+1);k++)
                 for(unsigned int j=0;j<(eleOrder+1);j++)
                   for(unsigned int i=0;i<(eleOrder+1);i++)
                   {
                       assert((sz%eleOrder)==0);
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



        void mesh2vtu(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames,const double **pointData)
        {


            if(!(pMesh->isActive())) return ;

            unsigned int rank=pMesh->getMPIRank();
            unsigned int npes=pMesh->getMPICommSize();

            unsigned int retval;


            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }

            const std::vector<ot::TreeNode> pElements=pMesh->getAllElements();
            const std::vector<unsigned int> e2N=pMesh->getE2NMapping();

            unsigned int dim=m_uiDim;
            DendroIntL num_vertices=pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement();
            unsigned int num_cells=pMesh->getNumLocalMeshElements();

            //int num_cells_elements = num_cells * NUM_CHILDREN + num_cells;
            unsigned int nPe=pMesh->getNumNodesPerElement();
            unsigned int eleOrder=pMesh->getElementOrder();
            unsigned int sz;

            std::vector<double> nodalVal;
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
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {


                            fprintf(fp, "          %d %d %d\n", (pElements[ele].getX() + i * (sz / eleOrder)),
                                    (pElements[ele].getY() + j * (sz / eleOrder)),
                                    (pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
#else

            unsigned int* coord_data = new unsigned int[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 0]=(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 1]=(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 2]=(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                  std::cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std::endl;
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
            DendroIntL * locidx_data =new DendroIntL [num_cells*NUM_CHILDREN];
            for (unsigned int i = 0 ; i < num_cells ; i++)
            {
                locidx_data[i*NUM_CHILDREN+0]=(DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+1]=(DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+2]=(DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+3]=(DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+4]=(DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                locidx_data[i*NUM_CHILDREN+5]=(DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+6]=(DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                locidx_data[i*NUM_CHILDREN+7]=(DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);

            }
            retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std::endl;
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
            DendroIntL * loc_offset=new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1]=NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std::endl;
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
            uint8_t  * loc_type=new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il]=VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std::endl;
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
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);

            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il]=rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_rank;


#endif
            fprintf (fp, "\n");
            fprintf(fp,"        </DataArray>\n");


#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pElements[il+pMesh->getElementLocalBegin()].getLevel());
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il]=pElements[il+pMesh->getElementLocalBegin()].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std::endl;
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
                double * nodalVal_binary=new double [num_cells*nPe];
#endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

#ifdef DENDRO_VTU_ASCII
                  fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
#else
                  fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
#endif


                  const double *tmp_var=pointData[pdata];
                  // note that size of tmp_var should be actual number of nodes in the mesh.
                  //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
#ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il);
                        for(unsigned int node=0;node<nPe;node++)
                            fprintf(fp,"%f ",nodalVal[node]);
                    }
                    fprintf(fp,"\n");
#else



                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il);
                        for(unsigned int node=0;node<nPe;node++)
                            nodalVal_binary[(il-pMesh->getElementLocalBegin())*nPe+node]=nodalVal[node];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_cells*nPe);

                    if (retval) {

                        std::cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std::endl;
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
            fp=NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp=fopen(fname,"w+");
                if(fp==NULL) {
                    std::cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std::endl;
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
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
#endif



#ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
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
                std::string vtuName=getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp=NULL;




            }


        }



        void mesh2vtuCoarse(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData)
        {

            if(!(pMesh->isActive())) return ;

            unsigned int rank=pMesh->getMPIRank();
            unsigned int npes=pMesh->getMPICommSize();

            int retval;

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }

            const std::vector<ot::TreeNode> pElements=pMesh->getAllElements();
            const std::vector<unsigned int> e2N=pMesh->getE2NMapping();

            unsigned int dim=m_uiDim;
            DendroIntL num_vertices=pMesh->getNumLocalMeshElements()*NUM_CHILDREN;
            unsigned int num_cells=pMesh->getNumLocalMeshElements();
            unsigned int nPe=pMesh->getNumNodesPerElement();
            unsigned int eleOrder=pMesh->getElementOrder();
            unsigned int sz;

            std::vector<double> nodalVal;
            nodalVal.resize(nPe);

            const unsigned int nx=2;
            const unsigned int ny=2;
            const unsigned int nz=2;

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
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k+=eleOrder)
                    for (unsigned int j = 0; j < (eleOrder + 1); j+=eleOrder)
                        for (unsigned int i = 0; i < (eleOrder + 1); i+=eleOrder) {


                            fprintf(fp, "          %d %d %d\n", (pElements[ele].getX() + i * (sz / eleOrder)),
                                    (pElements[ele].getY() + j * (sz / eleOrder)),
                                    (pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
#else

            unsigned int* coord_data = new unsigned int[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k+=eleOrder)
                    for (unsigned int j = 0; j < (eleOrder + 1); j+=eleOrder)
                        for (unsigned int i = 0; i < (eleOrder + 1); i+=eleOrder) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 0]=(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 1]=(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder))*m_uiDim + 2]=(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std::endl;
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


                DendroIntL * locidx_data =new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int i = 0 ; i < num_cells ; i++)
                {
                    locidx_data[i*NUM_CHILDREN+0]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+1]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+2]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+3]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0);
                    locidx_data[i*NUM_CHILDREN+4]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+5]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+6]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+7]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0);

                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std::cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std::endl;
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

            DendroIntL * loc_offset=new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1]=NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std::endl;
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
            uint8_t  * loc_type=new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il]=VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_type;

#endif

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il]=rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_rank;
#endif

            fprintf(fp,"        </DataArray>\n");


#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pElements[il+pMesh->getElementLocalBegin()].getLevel());
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il]=pElements[il+pMesh->getElementLocalBegin()].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std::endl;
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
                double * nodalVal_binary=new double [num_vertices];
#endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

#ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
#else
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
#endif


                    const double *tmp_var=pointData[pdata];
                    // note that size of tmp_var should be actual number of nodes in the mesh.
                    //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
#ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var, &(*(nodalVal.begin())), il);
                        for (unsigned int k = 0; k < (eleOrder + 1); k += eleOrder)
                            for (unsigned int j = 0; j < (eleOrder + 1); j += eleOrder)
                                for (unsigned int i = 0; i < (eleOrder + 1); i += eleOrder)
                                    fprintf(fp,"%f ",nodalVal[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]);
                    }
                    fprintf(fp,"\n");
#else



                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il);
                        for (unsigned int k = 0; k < (eleOrder + 1); k += eleOrder)
                            for (unsigned int j = 0; j < (eleOrder + 1); j += eleOrder)
                                for (unsigned int i = 0; i < (eleOrder + 1); i += eleOrder)
                                    nodalVal_binary[(il-pMesh->getElementLocalBegin())*NUM_CHILDREN+(k/eleOrder)*ny*nx+(j/eleOrder)*nx+(i/eleOrder)]=nodalVal[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_vertices);

                    if (retval) {

                        std::cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std::endl;
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
            fp=NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp=fopen(fname,"w+");
                if(fp==NULL) {
                    std::cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std::endl;
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
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
#endif



#ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
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
                std::string vtuName=getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp=NULL;




            }


        }


        void mesh2vtuFine(const ot::Mesh *pMesh, const char *fPrefix,unsigned int numFieldData,const char** filedDataNames,const double * filedData,unsigned int numPointdata, const char **pointDataNames, const double **pointData)
        {

            if(!(pMesh->isActive())) return ;

            unsigned int rank=pMesh->getMPIRank();
            unsigned int npes=pMesh->getMPICommSize();

            int retval;

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }

            const std::vector<ot::TreeNode> pElements=pMesh->getAllElements();
            const std::vector<unsigned int> e2N=pMesh->getE2NMapping();

            // note this works only for the 3D grids.
            const unsigned int dim=m_uiDim;
            const unsigned int nPe=pMesh->getNumNodesPerElement();
            const unsigned int eleOrder=pMesh->getElementOrder();
            const unsigned int ePe=eleOrder*eleOrder*eleOrder;
            const unsigned int num_cells=pMesh->getNumLocalMeshElements()*ePe;
            const DendroIntL num_vertices=pMesh->getNumLocalMeshElements()*nPe;
            unsigned int sz;

            std::vector<double> nodalVal;
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
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {


                            fprintf(fp, "          %d %d %d\n", (pElements[ele].getX() + i * (sz / eleOrder)),
                                    (pElements[ele].getY() + j * (sz / eleOrder)),
                                    (pElements[ele].getZ() + k * (sz / eleOrder)));

                        }
            }
#else

            unsigned int* coord_data = new unsigned int[num_vertices*m_uiDim];

            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++) {
                sz = 1u << (m_uiMaxDepth - pElements[ele].getLevel());
                assert((sz % eleOrder) == 0);
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 0]=(pElements[ele].getX()+i*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 1]=(pElements[ele].getY()+j*(sz/eleOrder));
                            coord_data[((ele-pMesh->getElementLocalBegin())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i)*m_uiDim + 2]=(pElements[ele].getZ()+k*(sz/eleOrder));
                        }
            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std::endl;
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

                DendroIntL * locidx_data =new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int ele = 0 ; ele < pMesh->getNumLocalMeshElements() ; ele++)
                {
                    for(unsigned int ek=0;ek<(eleOrder);ek++)
                        for(unsigned int ej=0;ej<(eleOrder);ej++)
                            for(unsigned int ei=0;ei<(eleOrder);ei++)
                            {

                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+0]=(DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+1]=(DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+2]=(DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+3]=(DendroIntL)(ele * nPe + ek*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+4]=(DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+ei);
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+5]=(DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+ej*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+6]=(DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+(ei+1));
                                locidx_data[(ele*ePe+ek*(eleOrder)*(eleOrder)+ej*(eleOrder)+ei)*NUM_CHILDREN+7]=(DendroIntL)(ele * nPe + (ek+1)*(eleOrder+1)*(eleOrder+1)+(ej+1)*(eleOrder+1)+ei);
                            }




                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std::cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std::endl;
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
            DendroIntL * loc_offset=new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1]=NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std::endl;
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
            uint8_t  * loc_type=new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il]=VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std::endl;
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
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il]=rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_rank;
#endif
            fprintf(fp,"\n");
            fprintf(fp,"        </DataArray>\n");


#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
               for(unsigned int w=0;w<ePe;w++)
                  fprintf(fp,"%d ",(pElements[il].getLevel() + (eleOrder>>1u)));
            }

            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
            unsigned int * loc_level=new unsigned int [num_cells];
            for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                for(unsigned int w=0;w<ePe;w++)
                    loc_level[(il-pMesh->getElementLocalBegin())*ePe+w]=pElements[il].getLevel()+(eleOrder>>1u);
            }

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std::endl;
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
                double * nodalVal_binary=new double [num_vertices];
#endif

                for(unsigned int pdata=0;pdata<numPointdata;pdata++)
                {

#ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_ASCII);
#else
                    fprintf(fp,"        <DataArray type=\"%s\" Name=\"%s\"" " format=\"%s\">\n",DENDRO_NODE_VAR_DOUBLE,pointDataNames[pdata],DENDRO_FORMAT_BINARY);
#endif


                    const double *tmp_var=pointData[pdata];
                    // note that size of tmp_var should be actual number of nodes in the mesh.
                    //std::cout<<"tmp_var: pdata: "<<pdata<<" var: "<<tmp_var[0]<<std::endl;
#ifdef DENDRO_VTU_ASCII
                    fprintf(fp,"         ");
                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il);
                        for(unsigned int node=0;node<nPe;node++)
                            fprintf(fp,"%f ",nodalVal[node]);
                    }
                    fprintf(fp,"\n");
#else



                    for (unsigned int il = pMesh->getElementLocalBegin(); il < pMesh->getElementLocalEnd(); ++il) {
                        pMesh->getElementNodalValues(tmp_var,&(*(nodalVal.begin())),il);
                        for(unsigned int node=0;node<nPe;node++)
                            nodalVal_binary[(il-pMesh->getElementLocalBegin())*nPe+node]=nodalVal[node];
                    }


                    retval = vtk_write_binary (fp, (char *) nodalVal_binary, sizeof (*nodalVal_binary)*num_vertices);

                    if (retval) {

                        std::cout<<rank<<": [VTU Error]: "<<"base64 encode point vars data failed"<<std::endl;
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
            fp=NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp=fopen(fname,"w+");
                if(fp==NULL) {
                    std::cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std::endl;
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
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
#endif



#ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
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
                std::string vtuName=getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp=NULL;




            }
        }


        void oct2vtu(const ot::TreeNode *pNodes,const unsigned int nSize,const char* fPrefix, MPI_Comm comm)
        {
            int rank,npes;
            MPI_Comm_rank(comm,&rank);
            MPI_Comm_size(comm,&npes);

            unsigned int retval;


            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }


            unsigned int dim=m_uiDim;
            DendroIntL num_vertices=nSize*NUM_CHILDREN;
            unsigned int num_cells=nSize;
            unsigned int sz;

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

            const unsigned int nx=2;
            const unsigned int ny=2;
            const unsigned int nz=2;


#ifdef DENDRO_VTU_ASCII
            for(unsigned int ele=0;ele<nSize;ele++) {
                sz = 1u << (m_uiMaxDepth - pNodes[ele].getLevel());
                assert((sz % eleOrder) == 0);
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
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 0]=(pNodes[ele].getX()+i*sz);
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 1]=(pNodes[ele].getY()+j*sz);
                         coord_data[((ele)*NUM_CHILDREN+k*ny*nx+j*nx+i)*m_uiDim + 2]=(pNodes[ele].getZ()+k*sz);
                     }


            }

            retval = vtk_write_binary (fp, (char *) coord_data, sizeof (*coord_data) * m_uiDim * num_vertices);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode point data failed"<<std::endl;
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


                DendroIntL * locidx_data =new DendroIntL [num_cells*NUM_CHILDREN];
                for (unsigned int i = 0 ; i < num_cells ; i++)
                {
                    locidx_data[i*NUM_CHILDREN+0]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+1]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+2]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+3]=(DendroIntL)(i * nx*ny*nz + 0*ny*nx+1*nx+0);
                    locidx_data[i*NUM_CHILDREN+4]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+0);
                    locidx_data[i*NUM_CHILDREN+5]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+0*nx+1);
                    locidx_data[i*NUM_CHILDREN+6]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+1);
                    locidx_data[i*NUM_CHILDREN+7]=(DendroIntL)(i * nx*ny*nz + 1*ny*nx+1*nx+0);

                }
                retval = vtk_write_binary (fp, (char *) locidx_data, sizeof (*locidx_data) * NUM_CHILDREN * num_cells);
                if (retval) {

                    std::cout<<rank<<": [VTU Error]: "<<"base64 encode connectivity data failed"<<std::endl;
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

            DendroIntL * loc_offset=new DendroIntL [num_cells];
            for (unsigned int il = 1; il <= (num_cells); il++)
                loc_offset[il-1]=NUM_CHILDREN*il;

            retval = vtk_write_binary (fp, (char *) loc_offset, sizeof (*loc_offset) *(num_cells));
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode offset data failed"<<std::endl;
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
            uint8_t  * loc_type=new uint8_t [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_type[il]=VTK_HEXAHEDRON;

            retval = vtk_write_binary (fp, (char *) loc_type, sizeof (*loc_type) *num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element type data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_type;

#endif

            fprintf(fp,"        </DataArray>\n");
            fprintf(fp,"      </Cells>\n");

            // writing cell data
            fprintf(fp,"      <CellData>\n");

#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf (fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",rank);
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_rank=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_rank[il]=rank;

            retval = vtk_write_binary (fp, (char *) loc_rank, sizeof (*loc_rank)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element rank data failed"<<std::endl;
                fclose(fp);
            }
            delete [] loc_rank;
#endif

            fprintf(fp,"        </DataArray>\n");


#ifdef DENDRO_VTU_ASCII
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            fprintf(fp,"         ");
            for (unsigned int il = 0; il < num_cells; ++il) {
                fprintf(fp,"%d ",pNodes[il].getLevel());
            }
            fprintf(fp,"\n");
#else
            fprintf(fp,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            unsigned int * loc_level=new unsigned int [num_cells];
            for (unsigned int il = 0; il < num_cells; ++il)
                loc_level[il]=pNodes[il].getLevel();

            retval = vtk_write_binary (fp, (char *) loc_level, sizeof (*loc_level)*num_cells);
            if (retval) {

                std::cout<<rank<<": [VTU Error]: "<<"base64 encode element level data failed"<<std::endl;
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
            fp=NULL;

            if(!rank) {
                sprintf(fname, "%s.pvtu", fPrefix);

                fp = fopen(fname, "w+");
                if (fp == NULL) {
                    std::cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std::endl;
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
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
#endif



#ifdef DENDRO_VTU_ASCII
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
#else
                fprintf(fp,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);
#endif

                fprintf(fp,"    </PCellData>\n");
                std::string vtuName=getFileName(std::string(fPrefix));
                for(unsigned int proc=0;proc<npes;proc++)
                {
                    fprintf(fp,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",vtuName.c_str(),proc,npes);

                }
                fprintf(fp,"</PUnstructuredGrid>\n");
                fprintf(fp,"</VTKFile>\n");


                fclose(fp);
                fp=NULL;


            }



        }


        void waveletsToVTU(ot::Mesh* mesh,const double ** zippedVec, const double **unzippedVec,const unsigned int * varIds,const unsigned int numVars,const char* fileName)
        {
    
            if(mesh->isActive())
            {
                ot::TreeNode blkNode;
                unsigned int sz[3];
                double dh[3];
                unsigned int bflag,offset;
                unsigned int regLev;
                unsigned int eIndex[3];
                double *  waveletR=new double[NUM_REFINE_WAVELET_COEF];
                //double *  waveletC=new double[NUM_COARSE_WAVELET_COEF]; 
                
                
                const unsigned int paddWidth=3;
                unsigned int eleIndexMin=0,eleIndexMax=0;
                const std::vector<ot::Block> localBlockList=mesh->getLocalBlockList();
                const std::vector<ot::TreeNode> allElements=mesh->getAllElements();
                const unsigned int eleOrder=mesh->getElementOrder();
                const unsigned int * e2n_cg=&(*(mesh->getE2NMapping().begin()));
                const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));
                const unsigned int nPe=mesh->getNumNodesPerElement();
                
                
                
                const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
                const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
                unsigned int nodeLookup_cg,nodeLookup_dg;
                unsigned int ownerID,ii_x,jj_y,kk_z;
                // allocate memory for wavelets
                double ** refineWavelets=new double*[numVars];
                for(unsigned int var=0;var<numVars;var++)
                    refineWavelets[var]=mesh->createVector<double>(100.0);
                
//                  for(unsigned int var=0;var<numVars;var++)
//                      for(unsigned int node=mesh->getNodeLocalBegin();node<mesh->getNodeLocalEnd();node++)
//                          refineWavelets[var][node]=zippedVec[varIds[var]][node];
               
                
                for(unsigned blk=0;blk<localBlockList.size();blk++)
                {
                    blkNode=localBlockList[blk].getBlockNode();
                    sz[0]=localBlockList[blk].getAllocationSzX();
                    sz[1]=localBlockList[blk].getAllocationSzY();
                    sz[2]=localBlockList[blk].getAllocationSzZ();

                    bflag=localBlockList[blk].getBlkNodeFlag();
                    offset=localBlockList[blk].getOffset();

                    regLev=localBlockList[blk].getRegularGridLev();
                    eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;

                    //if(bflag!=0) continue;
                        
                    
                    for(unsigned int ele=localBlockList[blk].getLocalElementBegin();ele<localBlockList[blk].getLocalElementEnd();ele++)
                    {

                        if((allElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)>=m_uiMaxDepth) continue;
                        
                        eIndex[0]=(allElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                        eIndex[1]=(allElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                        eIndex[2]=(allElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                        if((bflag &(1u<<OCT_DIR_LEFT)) && eIndex[0]==eleIndexMin)   continue;
                        if((bflag &(1u<<OCT_DIR_DOWN)) && eIndex[1]==eleIndexMin)   continue;
                        if((bflag &(1u<<OCT_DIR_BACK)) && eIndex[2]==eleIndexMin)   continue;

                        if((bflag &(1u<<OCT_DIR_RIGHT)) && eIndex[0]==eleIndexMax)  continue;
                        if((bflag &(1u<<OCT_DIR_UP)) && eIndex[1]==eleIndexMax)     continue;
                        if((bflag &(1u<<OCT_DIR_FRONT)) && eIndex[2]==eleIndexMax)  continue;

                        for(unsigned int var=0;var<numVars;var++)
                        {
                            computeRefineWavelets(unzippedVec[varIds[var]],offset,eleOrder,eIndex,paddWidth,sz,waveletR);
                            
//                             for(unsigned int k=0;k<(eleOrder+1);k++)
//                                for(unsigned int j=0;j<(eleOrder+1);j++)
//                                   for(unsigned int i=0;i<(eleOrder+1);i++)
//                                       if(!mesh->isNodeHanging(ele,i,j,k))refineWavelets[var][e2n[ele*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]]=normL2(waveletR,8);//normLInfty(waveletR,8);
                                  
                            /*refineWavelets[var][e2n_cg[ele*nPe+1*(eleOrder+1)*(eleOrder+1)+1*(eleOrder+1)+1]]=waveletR[0];
                            refineWavelets[var][e2n_cg[ele*nPe+1*(eleOrder+1)*(eleOrder+1)+1*(eleOrder+1)+3]]=waveletR[1];
                            
                            refineWavelets[var][e2n_cg[ele*nPe+1*(eleOrder+1)*(eleOrder+1)+3*(eleOrder+1)+1]]=waveletR[2];
                            refineWavelets[var][e2n_cg[ele*nPe+1*(eleOrder+1)*(eleOrder+1)+3*(eleOrder+1)+3]]=waveletR[3];
                            
                            refineWavelets[var][e2n_cg[ele*nPe+3*(eleOrder+1)*(eleOrder+1)+1*(eleOrder+1)+1]]=waveletR[4];
                            refineWavelets[var][e2n_cg[ele*nPe+3*(eleOrder+1)*(eleOrder+1)+1*(eleOrder+1)+3]]=waveletR[5];
                            
                            refineWavelets[var][e2n_cg[ele*nPe+3*(eleOrder+1)*(eleOrder+1)+3*(eleOrder+1)+1]]=waveletR[6];
                            refineWavelets[var][e2n_cg[ele*nPe+3*(eleOrder+1)*(eleOrder+1)+3*(eleOrder+1)+3]]=waveletR[7];*/
                            
                            //std::cout<<"waveletR[0]: "<<waveletR[0]<<std::endl;
                            
                        }

                    }
                    
                }
                
                
//                 for(unsigned int var=0;var<numVars;var++)
//                     for(unsigned int node=mesh->getNodeLocalBegin();node<mesh->getNodeLocalEnd();node++)
//                         refineWavelets[var][node]-=zippedVec[varIds[var]][node];
                
                for(unsigned int var=0;var<numVars;var++)
                    mesh->performGhostExchange(refineWavelets[var]);
                
                char** pDataNames=new char*[numVars];
                char *  pDataNamesAll[]={"U_ALPHA_WC","U_CHI_WC","U_K_WC","U_GT0_WC","U_GT1_WC","U_GT2_WC","U_BETA0_WC","U_BETA1_WC","U_BETA2_WC","U_B0_WC","U_B1_WC","U_B2_WC","U_SYMGT0_WC","U_SYMGT1_WC","U_SYMGT2_WC","U_SYMGT3_WC","U_SYMGT4_WC","U_SYMGT5_WC","U_SYMAT0_WC","U_SYMAT1_WC","U_SYMAT2_WC","U_SYMAT3_WC","U_SYMAT4_WC","U_SYMAT5_WC"};
                for(unsigned int var=0;var<numVars;var++)
                {
                    pDataNames[var]=pDataNamesAll[varIds[var]];
                }

                io::vtk::mesh2vtuFine(mesh,fileName,0,NULL,NULL,numVars,(const char **)pDataNames,(const double **)refineWavelets);
                
                for(unsigned int var=0;var<numVars;var++)
                {
                   delete [] refineWavelets[var]; 
                }
                   
                
                delete [] refineWavelets;
                delete [] waveletR;
                delete [] pDataNames;
            
                
                
            }
            
        }
        
        

    } // end of namespace vtk

} // end of namespace io
