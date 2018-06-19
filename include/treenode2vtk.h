//
// Created by milinda on 2/8/16.
//

#ifndef SFCSORTBENCH_TREENODE2VTK_H
#define SFCSORTBENCH_TREENODE2VTK_H
#include "TreeNode.h"
//#include "Point.h"
//#include "nodeAndValues.h"
#include <algorithm>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#define VTK_HEXAHEDRON 12

// convert one treenode to vtk file compatibel string
//std::string treeNodeTovtk(const TreeNode& T,int mpi_rank,int hindex);

//void treeNodesTovtk(std::vector<ot::TreeNode>& nodes,int mpi_rank,std::string vtk_file_name,bool hsorted=false);
template <typename T>
void treeNodesTovtk(std::vector<T>& nodes,int mpi_rank,std::string vtk_file_name,bool hsorted=false);

template <typename T>
void treeNodesTovtk(std::vector<T> &nodes, int mpi_rank, std::string vtk_file_name, bool hsorted) {

    //std::sort(nodes.begin(), nodes.end());
    if (!mpi_rank) std::cout << "writing mesh to VTK file: " << vtk_file_name << std::endl;
    std::ostringstream convert;

#ifdef  HILBERT_ORDERING
    convert << vtk_file_name << "_H_" << mpi_rank << ".vtk";
#else
    convert << vtk_file_name << "_M_" << mpi_rank << ".vtk";
#endif

    //convert << vtk_file_name << "_" << mpi_rank << ".vtk";
    vtk_file_name = convert.str();

    std::ofstream myfile;
    myfile.open(vtk_file_name.c_str());

    myfile << "# vtk DataFile Version 2.0" << std::endl;
    myfile << "DENDRO OCTREES" << std::endl;
    myfile << "ASCII" << std::endl;
    myfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int dim = 3;//nodes[0].getDim();

    int unit_points = 1 << dim;
    int num_verticies = nodes.size() * (unit_points);
    int num_cells = nodes.size();

    //std::cout<<"FileName :"<<vtk_file_name<<"\t DIM:"<<dim<<"\t number of points:"<<num_verticies<<"\t number of cells:"<<num_cells<<std::endl;


    myfile << "POINTS " << num_verticies << " float" << std::endl;

    if (dim == 2) {

    } else if (dim == 3) {
        unsigned int len;
        unsigned int xl, yl, zl;
        int num_data_field = 2;
        /* if (hsorted) {
             num_data_field++;
             std::sort(nodes.begin(), nodes.end());
         }*/

        int num_cells_elements = num_cells * unit_points + num_cells;

        for (int i = 0; i < nodes.size(); i++) {
            //std::cout<<nodes[i]<<std::endl;
            len = 1 << (nodes[i].getMaxDepth() - nodes[i].getLevel());
            xl = nodes[i].getX();
            yl = nodes[i].getY();
            zl = nodes[i].getZ();

            myfile << xl << " " << yl << " " << zl << std::endl;
            myfile << (xl + len) << " " << yl << " " << zl << std::endl;
            myfile << (xl + len) << " " << (yl + len) << " " << zl << std::endl;
            myfile << xl << " " << (yl + len) << " " << zl << std::endl;

            myfile << xl << " " << yl << " " << (zl + len) << std::endl;
            myfile << (xl + len) << " " << yl << " " << (zl + len) << std::endl;
            myfile << (xl + len) << " " << (yl + len) << " " << (zl + len) << std::endl;
            myfile << xl << " " << (yl + len) << " " << (zl + len) << std::endl;

        }

        myfile << "CELLS " << nodes.size() << " " << num_cells_elements << std::endl;

        for (int i = 0; i < num_cells; i++) {
            myfile << unit_points << " ";
            for (int j = 0; j < unit_points; j++) {
                myfile << (i * unit_points + j) << " ";
            }
            myfile << std::endl;
        }

        myfile << "CELL_TYPES " << num_cells << std::endl;
        for (int i = 0; i < num_cells; i++) {
            myfile << VTK_HEXAHEDRON << std::endl;
        }

        //myfile<<"CELL_DATA "<<num_cells<<std::endl;
        //myfile<<"POINT_DATA "<<(num_cells*unit_points)<<std::endl;



        myfile <<"CELL_DATA "<<num_cells<<std::endl;
        //myfile<<"POINT_DATA "<<(num_cells*unit_points)<<std::endl;



        //myfile << "FIELD OCTREE_DATA " << num_data_field << std::endl;

        //myfile << "cell_level 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS cell_level FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;

        for (int i = 0; i < nodes.size(); i++)
            myfile << nodes[i].getLevel() <<std::endl;

        myfile << std::endl;
        //myfile << std::endl;

        //myfile << "mpi_rank 1 " << num_cells << " int" << std::endl;
        myfile << "SCALARS mpi_rank FLOAT"<< std::endl;
        myfile << "LOOKUP_TABLE default"<<std::endl;
        for (int i = 0; i < nodes.size(); i++)
            myfile <<mpi_rank << std::endl;

        myfile << std::endl;


    }

    myfile.close();

}


#endif //SFCSORTBENCH_TREENODE2VTK_H
