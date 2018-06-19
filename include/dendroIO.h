//
// Created by milinda on 12/15/16.
//

/**
 * @author Milinda Fernando
 * @author Hari Sundar
 * @author Rahul Sampath
 *
 * @breif Constains all the IO related functionalities in dendro5
 *
 *
 * */

#ifndef SFCSORTBENCH_DENDROIO_H
#define SFCSORTBENCH_DENDROIO_H

#include "TreeNode.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "assert.h"


// Files to handle raw data.

namespace IO{

/**
    @author Rahul Sampath
    @brief Reads a list of points from a file
    @param filename the file name
    @param pts the points
    */
int readPtsFromFile(char* filename, std::vector<double>& pts);

/**
  @author Ilya Lashuk
  @brief Reads a list of points and corresponding values from a file
  @param filename the file name
  @param pts the points
  @param data the values
  */
int readDataPtsFromFile(char* filename, std::vector<double>& pts, std::vector<double>& ptVals);

/**
  @author Rahul Sampath
  @brief Writes a list of points to a file
  @param filename the file name
  @param pts the points
  */
int writePtsToFile(char* filename,  std::vector<double>& pts);

/**
  @author Ilya Lashuk
  @brief Writes a list of points and corresponding values from a file
  @param filename the file name
  @param pts the points
  @param data the values
  */
int writeDataPtsToFile(char* filename, std::vector<double>& pts, std::vector<double>& data);


/**
  @author Rahul Sampath
  @brief Writes a list of octants to a file
  @param filename the file name
  @param nodes the octants
  */
int writeNodesToFile (char* filename, const std::vector<ot::TreeNode> & nodes);

/**
  @author Rahul Sampath
  @brief Reads a list of octants from a file
  @param filename the file name
  @param nodes the octants
  */
int readNodesFromFile (char* filename,std::vector<ot::TreeNode > & nodes );



};



#endif //SFCSORTBENCH_DENDROIO_H
