//
// Created by milinda on 12/12/18.
//

/**
 * @brief contains functions to dump out raw data asci/binary. Can be useful for debugging purposes.
 * @author Milinda Fernando
 * */

#ifndef DENDRO_5_0_RAWIO_H
#define DENDRO_5_0_RAWIO_H

#include <iostream>
#include "mesh.h"
#include <fstream>
#include <TreeNode.h>


namespace io
{
    template <typename T>
    void varToRawData(const ot::Mesh* pMesh, const T** vars, unsigned int numVars,const char** varNames=NULL,const char* fPrefix="dendro_raw");

    template<typename T>
    void dump_array(const  T* const in, unsigned int sz, const char* fPrefix, MPI_Comm comm);



}// end of namespace io


#include "rawIO.tcc"


#endif //DENDRO_5_0_RAWIO_H
