//
// Created by milinda on 6/12/16.
//
/*
 * @author: Milinda Shayamal Fernando
 * School of Computing, University of Utah
 * @date: 6/12/2016
 *
 * This file contains code to test octree related utilities, like octree construction, balancing, implementations based on the treeSort aproach.
 *
 * */

#ifndef SFCSORTBENCH_TESTOCTUTILS_H
#define SFCSORTBENCH_TESTOCTUTILS_H

#include <iostream>
#include "treenode2vtk.h"
#include "genPts_par.h"
#include "TreeNode.h"
#include "mpi.h"
#include "dendro.h"
#include "sfcSort.h"
#include "treenode2vtk.h"
#include "testUtils.h"

#define ROOT_PROC 0



int readPtsFromFile(char* filename, std::vector<double>& pts);




#endif //SFCSORTBENCH_TESTOCTUTILS_H
