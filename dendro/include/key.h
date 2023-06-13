//
// Created by milinda on 6/28/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief contains treenode based search keys.
*/
//

#ifndef SFCSORTBENCH_KEY_H
#define SFCSORTBENCH_KEY_H

#include "TreeNode.h"
#include <iostream>
#include <vector>
#include <climits>
#include <memory>

namespace ot
{


    /**
     * @breif Key class to store some additional info to build nodelist.
     * */
    class Key : public ot::TreeNode
    {

    protected:
        std::vector<unsigned int >m_uiOwnerList;
        std::vector<unsigned int> m_uiStencilIndexWithDirection;  // first 3 bits contains the direction (OCT_DIR_LEFT,OCT_DIR_RIGHT etc) and the other 5 bit contains the stencil index.  Hence k_s max is 31.
        unsigned int m_uiSearchResult;


    public:
        Key()
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiLevel=0;
            m_uiSearchResult=0;

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();

        }

        Key(unsigned int px, unsigned int py, unsigned int pz, unsigned int plevel,unsigned int pDim,unsigned int pMaxDepth)/*: ot::TreeNode(px,py,pz,plevel,pDim,pMaxDepth)*/
        {
            m_uiX=px;
            m_uiY=py;
            m_uiZ=pz;
            m_uiLevel=plevel;
            m_uiSearchResult=0;

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();


        }

        Key(unsigned int pLevel, unsigned int pMaxDepth)
        {
            m_uiX=0;
            m_uiY=0;
            m_uiZ=0;
            m_uiSearchResult=0;
            m_uiLevel=pLevel;

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();


        }

        Key(const ot::TreeNode node)
        {
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiSearchResult=0;

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();


        }

        ~Key()
        {
            // need to delete the allocated variables

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();

        }


        inline void operator= ( const ot::TreeNode & node )

        {

            //std::cout<<"shared_ptr_ copy begin0 : "<<*this<< " = "<<node<<std::endl;
            //std::cout<<"shared_ptr_ copy begin1 : "<<*this<<" owner Size: "<<m_uiOwnerList.get()->size()<< " = "<<node<<" node ownerSize : "<<node.m_uiOwnerList.get()->size()<<std::endl;
            m_uiX=node.getX();
            m_uiY=node.getY();
            m_uiZ=node.getZ();
            m_uiLevel=node.getFlag();
            m_uiSearchResult=0;

            m_uiOwnerList.clear();
            m_uiStencilIndexWithDirection.clear();

        }


        inline void setSearchResult(unsigned int pIndex)
        {
            m_uiSearchResult=pIndex;
        }

        inline void addOwner(unsigned int ownerLocalID){
            m_uiOwnerList.push_back(ownerLocalID);
        }



        inline void addStencilIndexAndDirection(unsigned int index, unsigned int direction)
        {
           m_uiStencilIndexWithDirection.push_back(((index<<3) | direction));
        }


        inline  unsigned int getOwnerListSize(){return m_uiOwnerList.size();}
        inline  std::vector<unsigned int>*getOwnerList(){return &m_uiOwnerList; }
        inline  std::vector<unsigned int>*getStencilIndexDirectionList(){return &m_uiStencilIndexWithDirection; }

        inline unsigned int getSearchResult()const {return m_uiSearchResult;}


    };

}


#endif //SFCSORTBENCH_KEY_H
