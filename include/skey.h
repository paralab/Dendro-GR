//
// Created by milinda on 8/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief contains treenode based search keys.
*/
//

#ifndef SFCSORTBENCH_SEARCHKEY_H
#define SFCSORTBENCH_SEARCHKEY_H

#include "TreeNode.h"
#include <iostream>
#include <vector>
#include <climits>
#include <memory>

namespace ot
{


    /**
     * @breif SearchKey class to store some additional info to build nodelist.
     * */
    class SearchKey : public ot::TreeNode
    {

    protected:
        int m_uiOwner;
        int m_uiStencilIndexWithDirection;



    public:
        SearchKey()
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiLevel=0;
            m_uiOwner=-1;
            m_uiStencilIndexWithDirection=-1;


        }

        SearchKey(unsigned int px, unsigned int py, unsigned int pz, unsigned int plevel,unsigned int pDim,unsigned int pMaxDepth)/*: ot::TreeNode(px,py,pz,plevel,pDim,pMaxDepth)*/
        {
            m_uiX=px;
            m_uiY=py;
            m_uiZ=pz;
            m_uiLevel=plevel;
            m_uiOwner=-1;
            m_uiStencilIndexWithDirection=-1;

        }

        SearchKey(unsigned int pLevel, unsigned int pMaxDepth)
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiLevel=pLevel;
            m_uiOwner=-1;
            m_uiStencilIndexWithDirection=-1;

        }

        SearchKey(const ot::TreeNode node)
        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=-1;
            m_uiStencilIndexWithDirection=-1;

        }

        ~SearchKey()
        {
            // need to delete the allocated variables
        }

        inline void operator= (const SearchKey& node )
        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=node.m_uiOwner;
            m_uiStencilIndexWithDirection=node.m_uiStencilIndexWithDirection;
        }

        inline void operator= ( const ot::TreeNode & node )

        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiOwner=-1;
            m_uiStencilIndexWithDirection=-1;
        }

        inline void addOwner(unsigned int ownerLocalID){
           m_uiOwner=ownerLocalID;
        }


        inline void addStencilIndexAndDirection(unsigned int index, unsigned int direction)
        {
            m_uiStencilIndexWithDirection=((index<<3) | direction);
        }

        inline void addStencilIndexAndDirection(unsigned int direction)
        {
            m_uiStencilIndexWithDirection=direction;
        }


        inline  int getOwner(){return m_uiOwner;}
        inline  int getStencilIndexDirectionList(){return m_uiStencilIndexWithDirection;}


    };

}


#endif //SFCSORTBENCH_SEARCHKEY_H
