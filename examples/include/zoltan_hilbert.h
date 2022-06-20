
/*
 * @author: Majid Rasoeli
 * School of Computing, University of Utah

 *
 * Benchmark to access the zoltan implementation of Hilbert Curve.
 *
 * */


#ifndef SFCSORTBENCH_ZOLTAN_HILBERT_H
#define SFCSORTBENCH_ZOLTAN_HILBERT_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <ctime>
#include "zoltan.h"
#include <iostream>
#include "genPts_par.h"

//#define MESH_RANGE 1<<30    // points in the mesh will be generated randomly between 0 and MESH_RANGE

typedef struct
{
    ZOLTAN_ID_TYPE x;
    ZOLTAN_ID_TYPE y;
    ZOLTAN_ID_TYPE z;
    ZOLTAN_ID_TYPE globalID;

} Node_Type;


/* Structure to hold mesh data */

typedef struct{
    ZOLTAN_ID_TYPE numGlobalPoints;
    ZOLTAN_ID_TYPE numLocalPts;
    ZOLTAN_ID_TYPE *globalIds;
    ZOLTAN_ID_TYPE *localIds;
    ZOLTAN_ID_TYPE *x;
    ZOLTAN_ID_TYPE *y;
    ZOLTAN_ID_TYPE *z;
} MESH_DATA;






void ZoltanHibertLBandSort(DendroIntL grainSz,unsigned int dim ,unsigned int maxDepth,double tol, std::vector<double>& stats,unsigned int distribution,MPI_Comm comm);

int get_number_of_objects(void *data, int *ierr);
void get_object_list(void *data, int sizeGID, int sizeLID,
                            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int wgt_dim, float *obj_wgts, int *ierr);
int get_num_geometry(void *data, int *ierr);
void get_geometry_list(void *data, int sizeGID, int sizeLID,
                              int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr);


void user_sizes_node(void *data,int num_gid_entries,int num_lid_entries, int nids, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int* nbytes, int *ierr);
void user_pack_node(void *data,int num_gid_entries,int num_lid_entries, int nids, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int* dest, int* size, int* index, char *buf, int *ierr);
void user_unpack_node(void *data,int gidSize,int num_ids, ZOLTAN_ID_PTR global_id, int * size,int * idx , char *buf, int *ierr);




namespace par {

    //Forward Declaration
    template <typename T>
    class Mpi_datatype;

/**
@author Rahul Sampath, rahul.sampath@gmail.com
@brief A template specialization of the abstract class "Mpi_datatype" for communicating messages of type "ot::TreeNode".
*/
    template <>
    class Mpi_datatype< Node_Type > {
/*

        static void Node_MAX_LEVEL(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] =
                        ( ( (first.getLevel()) > (second.getLevel()) )? first : second );
            }//end for
        }//end function

        static void Node_MAX(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first > second )? first : second );
            }//end for
        }//end function

        static void Node_MIN(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first < second )? first : second );
            }//end for
        }//end function
*/

        /* static void Node_NCA(void *in, void *inout, int* len, MPI_Datatype * dptr) {
             for(int i = 0; i < (*len); i++) {
                 ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                 ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                 if( first != second ) {
                     (static_cast<ot::TreeNode*>(inout))[i] = ot::getNCA(first, second);
                 }//end if
             }//end for
         }//end function*/

    public:
        /**
          @brief User defined MPI_Operation that sets second[i] to first[i] if first[i] is at a greater level than second[i].
          @remark first and second are 2 arrays of type TreeNode.
        **/
        /*static MPI_Op MAX_LEVEL(){
            static bool first = true;
            static MPI_Op maxLev;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX_LEVEL ,true ,&maxLev);
            }
            return maxLev;
        }*/

        /**
          @brief User defined MPI_Operation that computes: second[i] = Max(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode.
          "MAX" is a macro
        **/
        /*static MPI_Op _MAX() {
            static bool         first = true;
            static MPI_Op max;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX ,true ,&max);
            }
            return max;
        }*/

        /**
          @brief User defined MPI_Operation that computes: second[i] = Min(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode.
          "MIN" is a macro
        **/
        /*static MPI_Op _MIN() {
            static bool         first = true;
            static MPI_Op min;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MIN ,true ,&min);
            }
            return min;
        }*/

        /**
          @brief User defined MPI_Operation that computes: second[i] = NCA(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode and
          NCA returns the nearest common ancestor of its 2 arguments.
        **/
        /*static MPI_Op NCA() {
            static bool         first = true;
            static MPI_Op nca;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_NCA ,true ,&nca);
            }
            return nca;
        }*/

        /**
         @return The MPI_Datatype corresponding to the datatype "ot::TreeNode".
       */
        static MPI_Datatype value()
        {
            static bool         first = true;
            static MPI_Datatype datatype;

            if (first)
            {
                first = false;
                MPI_Type_contiguous(sizeof(Node_Type), MPI_BYTE, &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }

    };

}//end namespace par



#endif //SFCSORTBENCH_ZOLTAN_HILBERT_H
