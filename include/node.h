//
// Created by milinda on 9/11/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief class derived from the TreeNode for scatter map computations.
*/
//

#ifndef SFCSORTBENCH_NODE_H
#define SFCSORTBENCH_NODE_H

#include "TreeNode.h"
#include <iostream>
#include <vector>
#include <climits>
#include <memory>

namespace ot {


    /**
     * @breif SearchKey class to store some additional info to build nodelist.
     * */
    class Node : public ot::TreeNode {

    protected:
        int m_uiOwner; // set to -1 by default.

    public:
        Node()
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiLevel=0;
            m_uiOwner=-1;

        }

        Node(unsigned int px, unsigned int py, unsigned int pz, unsigned int plevel,unsigned int pDim,unsigned int pMaxDepth)/*: ot::TreeNode(px,py,pz,plevel,pDim,pMaxDepth)*/
        {
            m_uiX=px;
            m_uiY=py;
            m_uiZ=pz;
            m_uiLevel=plevel;
            m_uiOwner=-1;

        }

        Node(unsigned int pLevel, unsigned int pMaxDepth)
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiLevel=pLevel;
            m_uiOwner=-1;
        }

        Node(const ot::TreeNode node)
        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=-1;

        }

        ~Node()
        {
            // need to delete the allocated variables
        }

        inline void setOwner(unsigned int ownerID) {m_uiOwner=ownerID;}

        inline void operator= (const Node& node )
        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=node.m_uiOwner;

        }

        inline void operator= ( const ot::TreeNode & node )

        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=-1;

        }
        inline  int getOwner(){return m_uiOwner;}

    };

}


namespace par {

    //Forward Declaration
    template <typename T>
    class Mpi_datatype;

/**
@author Rahul Sampath, rahul.sampath@gmail.com
@brief A template specialization of the abstract class "Mpi_datatype" for communicating messages of type "ot::Node".
*/
    template <>
    class Mpi_datatype< ot::Node > {

//        static void Node_MAX_LEVEL(void *in, void *inout, int* len, MPI_Datatype * dptr) {
//            for(int i = 0; i < (*len); i++) {
//                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
//                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
//                (static_cast<ot::TreeNode*>(inout))[i] =
//                        ( ( (first.getLevel()) > (second.getLevel()) )? first : second );
//            }//end for
//        }//end function
//
//        static void Node_MAX(void *in, void *inout, int* len, MPI_Datatype * dptr) {
//            for(int i = 0; i < (*len); i++) {
//                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
//                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
//                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first > second )? first : second );
//            }//end for
//        }//end function
//
//        static void Node_MIN(void *in, void *inout, int* len, MPI_Datatype * dptr) {
//            for(int i = 0; i < (*len); i++) {
//                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
//                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
//                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first < second )? first : second );
//            }//end for
//        }//end function
//
//        /* static void Node_NCA(void *in, void *inout, int* len, MPI_Datatype * dptr) {
//             for(int i = 0; i < (*len); i++) {
//                 ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
//                 ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
//                 if( first != second ) {
//                     (static_cast<ot::TreeNode*>(inout))[i] = ot::getNCA(first, second);
//                 }//end if
//             }//end for
//         }//end function*/

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
//        static MPI_Op _MAX() {
//            static bool         first = true;
//            static MPI_Op max;
//            if (first) {
//                first = false;
//                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX ,true ,&max);
//            }
//            return max;
//        }

        /**
          @brief User defined MPI_Operation that computes: second[i] = Min(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode.
          "MIN" is a macro
        **/
//        static MPI_Op _MIN() {
//            static bool         first = true;
//            static MPI_Op min;
//            if (first) {
//                first = false;
//                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MIN ,true ,&min);
//            }
//            return min;
//        }

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
                MPI_Type_contiguous(sizeof(ot::Node), MPI_BYTE, &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }

    };

}//end namespace par


#endif //SFCSORTBENCH_NODE_H
