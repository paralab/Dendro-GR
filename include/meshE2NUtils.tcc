//
// Created by milinda on 9/22/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains the Element to nodal map generation helper functions for
 * octree based mesh generation. (these functions are transformed form the mesh.h file. )
*/
//

namespace ot
{
    inline void Mesh::OCT_DIR_LEFT_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex)
    {
        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;
        // -- EDGE OF Element e , (OCT_DIR_LEFT, OCT_DIR_UP)
        edgeOwner=parent;

        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_LEFT,OCT_DIR_UP,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);

        }
        else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_UP,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // -- EDGE OF Element e , (OCT_DIR_LEFT,OCT_DIR_DOWN)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}


        edgeNodeIndex(child,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeChildIndex,true);


        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);
        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }

        //-----------------------------------------------------


        // EDGE of Element e (OCT_DIR_LEFT, OCT_DIR_FRONT)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);
        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);

            }


            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }




        // EDGE of Element e (OCT_DIR_LEFT, OCT_DIR_BACK)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_LEFT,OCT_DIR_BACK,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);
        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeChildIndex,true);

            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);

            }

            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


    }

    inline void Mesh::OCT_DIR_RIGHT_INTERNAL_EDGE_MAP(unsigned int child, unsigned int parent, bool parentChildLevEqual, std::vector<unsigned int> &edgeChildIndex, std::vector<unsigned int> &edgeOwnerIndex)
    {
        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;
        // Edge of element e (OCT_DIR_RIGHT,OCT_DIR_UP)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_RIGHT,OCT_DIR_UP,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);


        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_UP,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // Edge of element e (OCT_DIR_RIGHT,OCT_DIR_DOWN)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];



        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_UP,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_UP,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


        // Edge of element e (OCT_DIR_RIGHT,OCT_DIR_FRONT)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


        // Edge of element e (OCT_DIR_RIGHT,OCT_DIR_BACK)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_BACK,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];



        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_LEFT,OCT_DIR_BACK,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


    }


    inline void Mesh::OCT_DIR_DOWN_INTERNAL_EDGE_MAP(unsigned int child, unsigned int parent, bool parentChildLevEqual, std::vector<unsigned int> &edgeChildIndex, std::vector<unsigned int> &edgeOwnerIndex)
    {

        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;

        // Edge of element e (OCT_DIR_DOWN,OCT_DIR_RIGHT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_RIGHT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


        // Edge of element e (OCT_DIR_DOWN,OCT_DIR LEFT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_LEFT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // Edge of element e (OCT_DIR_DOWN,OCT_DIR_FRONT)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];



        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_FRONT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // Edge of element e (OCT_DIR_DOWN,OCT_DIR_BACK)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        // }

        edgeNodeIndex(child,OCT_DIR_DOWN,OCT_DIR_BACK,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);


        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_UP,OCT_DIR_BACK,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }

    }

    inline void Mesh::OCT_DIR_UP_INTERNAL_EDGE_MAP(unsigned int child, unsigned int parent, bool parentChildLevEqual, std::vector<unsigned int> &edgeChildIndex, std::vector<unsigned int> &edgeOwnerIndex)
    {

        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;


        // Edge of element e (OCT_DIR_UP,OCT_DIR_LEFT)


        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_UP,OCT_DIR_LEFT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_RIGHT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }

        // Edge of Element e (OCT_DIR_UP, OCT_DIR_RIGHT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_UP,OCT_DIR_RIGHT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_RIGHT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_LEFT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_LEFT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // Edge of element e (OCT_DIR_UP,OCT_DIR_FRONT)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_UP,OCT_DIR_FRONT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_BACK,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


        // Edge of element e (OCT_DIR_UP,OCT_DIR_BACK)
        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        // {
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_UP,OCT_DIR_BACK,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_BACK,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_DOWN,OCT_DIR_BACK,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_UP,OCT_DIR_FRONT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


    }

    inline void Mesh::OCT_DIR_BACK_INTERNAL_EDGE_MAP(unsigned int child, unsigned int parent, bool parentChildLevEqual, std::vector<unsigned int> &edgeChildIndex, std::vector<unsigned int> &edgeOwnerIndex)
    {

        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;

        // edge of the element e (OCT_DIR_BACK, OCT_DIR_LEFT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_BACK,OCT_DIR_LEFT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }

        // edge of the element e (OCT_DIR_BACK, OCT_DIR_RIGHT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


        // edge of the element e (OCT_DIR_BACK, OCT_DIR_DOWN)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_BACK,OCT_DIR_DOWN,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }

        // edge of the element e (OCT_DIR_BACK, OCT_DIR_UP)


        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_BACK,OCT_DIR_UP,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_FRONT,OCT_DIR_UP,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



    }

    inline void Mesh::OCT_DIR_FRONT_INTERNAL_EDGE_MAP(unsigned int child, unsigned int parent, bool parentChildLevEqual, std::vector<unsigned int> &edgeChildIndex, std::vector<unsigned int> &edgeOwnerIndex)
    {

        unsigned int edgeOwner;
        unsigned int lookUp1,lookUp2;

        // edge of the element e (OCT_DIR_FRONT, OCT_DIR_LEFT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_LEFT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // edge of the element e (OCT_DIR_FRONT, OCT_DIR_RIGHT)

        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        // if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_FRONT,OCT_DIR_RIGHT,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];




        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_RIGHT,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_LEFT,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_LEFT,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }





        // edge of the element e (OCT_DIR_FRONT, OCT_DIR_DOWN)


        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        // if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);

        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];



        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_DOWN,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_UP,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }



        // edge of the element e (OCT_DIR_FRONT, OCT_DIR_UP)


        edgeOwner=parent;
        lookUp1=m_uiE2EMapping[parent*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp1!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp1].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp1].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp1<edgeOwner)))
                edgeOwner=lookUp1;
        }

        //if(parentChildLevEqual)
        //{
        lookUp2=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        if(lookUp2!=LOOK_UP_TABLE_DEFAULT)
        {
            if((m_uiAllElements[lookUp2].getLevel()< m_uiAllElements[edgeOwner].getLevel()) || ((m_uiAllElements[lookUp2].getLevel()==m_uiAllElements[edgeOwner].getLevel()) && (lookUp2<edgeOwner)))
                edgeOwner=lookUp2;
        }
        //}

        edgeNodeIndex(child,OCT_DIR_FRONT,OCT_DIR_UP,edgeChildIndex,true);

        if(edgeOwner==parent)
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_UP,edgeOwnerIndex,true);


        }else if(edgeOwner==lookUp1)
        {
            assert(m_uiAllElements[lookUp1].getLevel()<=m_uiAllElements[child].getLevel());
            if(lookUp2==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);
            else
                edgeNodeIndex(lookUp1,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);


        }else
        {
            assert(edgeOwner==lookUp2);
            edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);
        }

        assert(edgeChildIndex.size()==edgeOwnerIndex.size());
        for(unsigned int index=0;index<edgeChildIndex.size();index++)
            m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];


        if((parentChildLevEqual) && (edgeOwner!=parent) )
        {
            edgeNodeIndex(parent,OCT_DIR_BACK,OCT_DIR_UP,edgeChildIndex,true);
            if(edgeOwner==lookUp2)
            {
                assert(m_uiAllElements[lookUp2].getLevel()<=m_uiAllElements[parent].getLevel());
                if(lookUp1==LOOK_UP_TABLE_DEFAULT || m_uiAllElements[lookUp1]!=m_uiAllElements[lookUp2])
                    edgeNodeIndex(lookUp2,OCT_DIR_FRONT,OCT_DIR_DOWN,edgeOwnerIndex,true);
                else
                    edgeNodeIndex(lookUp2,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);
            }else
            {
                assert(edgeOwner==lookUp1);
                edgeNodeIndex(lookUp1,OCT_DIR_BACK,OCT_DIR_DOWN,edgeOwnerIndex,true);

            }
            assert(edgeChildIndex.size()==edgeOwnerIndex.size());
            for(unsigned int index=0;index<edgeChildIndex.size();index++)
                m_uiE2NMapping_CG[edgeChildIndex[index]]=m_uiE2NMapping_CG[edgeOwnerIndex[index]];

        }


    }


    inline void Mesh::CORNER_NODE_MAP(unsigned int child)
    {


        // DEBUG HELPER CODE TO figure out Nodal mismatches. Be really carefull when changing this function :)
        /*if(m_uiActiveRank==1) {
            std::vector<ot::TreeNode> cusEleCheck;
            cusEleCheck.push_back(m_uiAllElements[child]);
            unsigned int ownerID, ii_x, jj_y, kk_z, x, y, z, sz;
            dg2eijk(cornerOwnerIndex,ownerID,ii_x,jj_y,kk_z);
            x = m_uiAllElements[ownerID].getX();
            y = m_uiAllElements[ownerID].getY();
            z = m_uiAllElements[ownerID].getZ();
            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
            cusEleCheck.push_back(
                    ot::TreeNode((x + ii_x * sz / m_uiElementOrder), (y + jj_y * sz / m_uiElementOrder),
                                 (z + kk_z * sz / m_uiElementOrder), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            treeNodesTovtk(cusEleCheck,child,"cusEleCheck");
        }*/


        unsigned int cornerOwner;
        unsigned int cornerOwnerIndex;
        unsigned int cornerChildIndex;
        unsigned int numLookUp=6;
        unsigned int lookUp[numLookUp];

        // MORTON INDEX 0 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,0,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_LEFT,OCT_DIR_BACK,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_DOWN,OCT_DIR_LEFT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_BACK,OCT_DIR_DOWN,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }

        assert(cornerOwner<m_uiAllElements.size());

        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],1,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],2,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],4,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],5,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],3,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],6,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];


        // MORTON INDEX 1 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,1,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_RIGHT,OCT_DIR_BACK,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_DOWN,OCT_DIR_RIGHT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_BACK,OCT_DIR_DOWN,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],0,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],3,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],5,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],4,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],2,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],7,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];


        // MORTON INDEX 2 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,2,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_LEFT,OCT_DIR_BACK,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_UP,OCT_DIR_LEFT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_BACK,OCT_DIR_UP,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],3,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],0,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],6,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],7,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],1,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],4,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];



        // MORTON INDEX 3 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,3,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_BACK];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_RIGHT,OCT_DIR_BACK,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_UP,OCT_DIR_RIGHT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_BACK,OCT_DIR_UP,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],2,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],1,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],7,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],6,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],0,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],5,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];


        // MORTON INDEX 4 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,4,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_LEFT,OCT_DIR_FRONT,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_DOWN,OCT_DIR_LEFT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_FRONT,OCT_DIR_DOWN,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],5,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],6,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],0,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],1,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],7,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],2,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];


        // MORTON INDEX 5 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,5,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_DOWN];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_RIGHT,OCT_DIR_FRONT,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_DOWN,OCT_DIR_RIGHT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_FRONT,OCT_DIR_DOWN,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],4,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],7,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],1,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],0,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],6,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],3,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];


        // MORTON INDEX 6 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,6,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_LEFT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_LEFT,OCT_DIR_FRONT,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_UP,OCT_DIR_LEFT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_FRONT,OCT_DIR_UP,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],7,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {

            cornerNodeIndex(lookUp[1],4,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],2,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],3,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],5,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],0,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];



        // MORTON INDEX 7 NODE.
        cornerOwner=child;
        cornerNodeIndex(child,7,cornerChildIndex);

        lookUp[0]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_RIGHT];
        lookUp[1]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_UP];
        lookUp[2]=m_uiE2EMapping[child*m_uiNumDirections+OCT_DIR_FRONT];

        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_RIGHT,OCT_DIR_FRONT,lookUp[3]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_UP,OCT_DIR_RIGHT,lookUp[4]);
        OCT_DIR_DIAGONAL_E2E(child,OCT_DIR_FRONT,OCT_DIR_UP,lookUp[5]);


        for(unsigned int w=0;w<numLookUp;w++)
        {
            if( (lookUp[w]!=LOOK_UP_TABLE_DEFAULT ) && ((m_uiAllElements[lookUp[w]].getLevel()<m_uiAllElements[cornerOwner].getLevel()) || ((m_uiAllElements[lookUp[w]].getLevel()==m_uiAllElements[cornerOwner].getLevel()) && (lookUp[w]<cornerOwner))))
                cornerOwner=lookUp[w];
        }


        if(cornerOwner==lookUp[0])
        {
            cornerNodeIndex(lookUp[0],6,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[1])
        {
            cornerNodeIndex(lookUp[1],5,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[2])
        {
            cornerNodeIndex(lookUp[2],3,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[3])
        {
            cornerNodeIndex(lookUp[3],2,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[4])
        {
            cornerNodeIndex(lookUp[4],4,cornerOwnerIndex);

        }else if(cornerOwner==lookUp[5])
        {
            cornerNodeIndex(lookUp[5],1,cornerOwnerIndex);
        }

        if(cornerOwner!=child) m_uiE2NMapping_CG[cornerChildIndex]=m_uiE2NMapping_CG[cornerOwnerIndex];



    }




    inline void Mesh::faceNodesIndex(unsigned int elementID, unsigned int face, std::vector<unsigned int>& index, bool isInternal) const {


        unsigned int begin=0;
        unsigned int end=(m_uiElementOrder+1);
        unsigned int numPts1D;

        if(isInternal) {
            begin = 1;
            end = (m_uiElementOrder);
        }
        numPts1D=(end-begin);
        if(numPts1D==0)
        {
            index.clear();
            return;
        }

        index.resize(numPts1D*numPts1D);


        if(face ==OCT_DIR_LEFT)
        {
            for(unsigned int k=begin;k<end;k++)
            {
                for(unsigned int j=begin;j<end;j++)
                {
                    index[(k-begin)*(numPts1D)+(j-begin)]=elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0;
                }
            }

        }else if(face==OCT_DIR_RIGHT)
        {
            for(unsigned int k=begin;k<end;k++)
            {
                for(unsigned int j=begin;j<end;j++)
                {
                    index[(k-begin)*(numPts1D)+(j-begin)]=elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder;
                }
            }


        }else if(face==OCT_DIR_DOWN)
        {

            for(unsigned int k=begin;k<end;k++)
            {
                for(unsigned int i=begin;i<end;i++)
                {
                    index[(k-begin)*(numPts1D)+(i-begin)]=elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i;
                }
            }

        }else if(face==OCT_DIR_UP)
        {

            for(unsigned int k=begin;k<end;k++)
            {
                for(unsigned int i=begin;i<end;i++)
                {
                    index[(k-begin)*(numPts1D)+(i-begin)]=elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i;
                }
            }

        }else if(face==OCT_DIR_FRONT)
        {

            for(unsigned int j=begin;j<end;j++)
            {
                for(unsigned int i=begin;i<end;i++)
                {
                    index[(j-begin)*(numPts1D)+(i-begin)]=elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i;
                }
            }

        }else if(face==OCT_DIR_BACK)
        {

            for(unsigned int j=begin;j<end;j++)
            {
                for(unsigned int i=begin;i<end;i++)
                {
                    index[(j-begin)*(numPts1D)+(i-begin)]=elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i;
                }
            }

        }


    }


    inline void Mesh::edgeNodeIndex(unsigned int elementID,unsigned int face1,unsigned int face2, std::vector<unsigned int >& index, bool isInternal) const
    {
        unsigned int begin=0;
        unsigned int end=(m_uiElementOrder+1);
        unsigned int numPts1D;



        if(isInternal) {
            begin = 1;
            end = (m_uiElementOrder);
        }
        numPts1D=(end-begin);
        if(numPts1D==0)
        {
            index.clear();
            return;
        }

        index.resize(numPts1D);
        assert(face1!=face2);
        //std::cout<<"begin: "<<begin<<" end: "<<end<<std::endl;

        if(((face1==OCT_DIR_LEFT)|| (face1==OCT_DIR_UP)) && ( (face2==OCT_DIR_UP) || (face2==OCT_DIR_LEFT) ))
        {
            for(unsigned int k=begin;k<end;k++)
                index[(k-begin)] = elementID*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + m_uiElementOrder*(m_uiElementOrder+1) + 0;

        }else if(((face1==OCT_DIR_LEFT)|| (face1==OCT_DIR_DOWN)) && ( (face2==OCT_DIR_DOWN) || (face2==OCT_DIR_LEFT) ))
        {
            for(unsigned int k=begin;k<end;k++)
                index[(k-begin)] = elementID*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + 0*(m_uiElementOrder+1) + 0;

        }else if(((face1==OCT_DIR_LEFT)|| (face1==OCT_DIR_FRONT)) && ( (face2==OCT_DIR_FRONT) || (face2==OCT_DIR_LEFT) ))
        {
            for(unsigned int j=begin;j<end;j++)
                index[(j-begin)] = elementID*m_uiNpE + m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + 0;

        }else if(((face1==OCT_DIR_LEFT)|| (face1==OCT_DIR_BACK)) && ( (face2==OCT_DIR_BACK) || (face2==OCT_DIR_LEFT) ))
        {
            for(unsigned int j=begin;j<end;j++)
                index[(j-begin)] = elementID*m_uiNpE + 0*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + 0;

        }else if(((face1==OCT_DIR_RIGHT)|| (face1==OCT_DIR_UP)) && ( (face2==OCT_DIR_UP) || (face2==OCT_DIR_RIGHT) )) {

            for(unsigned int k=begin;k<end;k++)
                index[(k-begin)] = elementID*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + m_uiElementOrder*(m_uiElementOrder+1) + m_uiElementOrder;

        }else if(((face1==OCT_DIR_RIGHT)|| (face1==OCT_DIR_DOWN)) && ( (face2==OCT_DIR_DOWN) || (face2==OCT_DIR_RIGHT) ))
        {
            for(unsigned int k=begin;k<end;k++)
                index[(k-begin)] = elementID*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + 0*(m_uiElementOrder+1) + m_uiElementOrder;

        }else if(((face1==OCT_DIR_RIGHT)|| (face1==OCT_DIR_FRONT)) && ( (face2==OCT_DIR_FRONT) || (face2==OCT_DIR_RIGHT) ))
        {
            for(unsigned int j=begin;j<end;j++)
                index[(j-begin)] = elementID*m_uiNpE + m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + m_uiElementOrder;

        }else if(((face1==OCT_DIR_RIGHT)|| (face1==OCT_DIR_BACK)) && ( (face2==OCT_DIR_BACK) || (face2==OCT_DIR_RIGHT) ))
        {
            for(unsigned int j=begin;j<end;j++)
                index[(j-begin)] = elementID*m_uiNpE + 0*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + m_uiElementOrder;

        }else if(((face1==OCT_DIR_UP)|| (face1==OCT_DIR_FRONT)) && ( (face2==OCT_DIR_FRONT) || (face2==OCT_DIR_UP) ))
        {
            for(unsigned int i=begin;i<end;i++)
                index[(i-begin)] = elementID*m_uiNpE + m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1) + m_uiElementOrder*(m_uiElementOrder+1)+ i;

        }else if(((face1==OCT_DIR_UP)|| (face1==OCT_DIR_BACK)) && ( (face2==OCT_DIR_BACK) || (face2==OCT_DIR_UP) ))
        {
            for(unsigned int i=begin;i<end;i++)
                index[(i-begin)] = elementID*m_uiNpE +0*(m_uiElementOrder+1)*(m_uiElementOrder+1) + m_uiElementOrder*(m_uiElementOrder+1) + i;

        }else if(((face1==OCT_DIR_DOWN)|| (face1==OCT_DIR_FRONT)) && ( (face2==OCT_DIR_FRONT) || (face2==OCT_DIR_DOWN) ))
        {
            for(unsigned int i=begin;i<end;i++)
                index[(i-begin)] = elementID*m_uiNpE + m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1) + 0*(m_uiElementOrder+1) + i;

        }else if(((face1==OCT_DIR_DOWN)|| (face1==OCT_DIR_BACK)) && ( (face2==OCT_DIR_BACK) || (face2==OCT_DIR_DOWN) ))
        {
            for(unsigned int i=begin;i<end;i++)
                index[(i-begin)] = elementID*m_uiNpE + 0*(m_uiElementOrder+1)*(m_uiElementOrder+1) + 0*(m_uiElementOrder+1) + i;
        }





    }

    inline void Mesh::cornerNodeIndex(unsigned int elementID,unsigned int mortonIndex, unsigned int &index) const
    {
        if(mortonIndex==0)
        {
            index=elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0;
        }else if(mortonIndex==1)
        {
            index=elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder;
        }else if(mortonIndex==2)
        {
            index=elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0;

        }else if(mortonIndex==3)
        {
            index=elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder;

        }else if(mortonIndex==4)
        {
            index=elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0;

        }else if(mortonIndex==5)
        {
            index=elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder;

        }else if(mortonIndex==6)
        {
            index=elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0;

        }else if(mortonIndex==7)
        {
            index=elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder;
        }


    }


    inline void Mesh::elementNodeIndex(unsigned int elementID,std::vector<unsigned int > & index, bool isInternal) const
    {
        unsigned int begin=0;
        unsigned int end=(m_uiElementOrder+1);
        unsigned int numPts1D;

        if(isInternal) {
            begin = 1;
            end = (m_uiElementOrder);
        }
        numPts1D=(end-begin);
        if(numPts1D==0)
        {
            index.clear();
            return;
        }

        index.resize(numPts1D*numPts1D*numPts1D);

        for(unsigned int k=begin;k<end;k++)
        {
            for(unsigned int j=begin;j<end;j++)
            {
                for(unsigned int i=begin;i<end;i++)
                {
                    index[(k-begin)*numPts1D*numPts1D+(j-begin)*numPts1D+(i-begin)]=elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + i;
                }
            }
        }


    }



    inline void Mesh::OCT_DIR_DIAGONAL_E2E(unsigned int elementID, unsigned int face1,unsigned int face2,unsigned int& lookUp) const
    {

        assert(elementID<m_uiAllElements.size());
        unsigned int myLev=m_uiAllElements[elementID].getLevel();
        /* unsigned int lookUpF1;
         unsigned int lookUpF2;*/

        lookUp=m_uiE2EMapping[elementID*m_uiNumDirections+face1];

        /*   if(m_uiActiveRank==0 && elementID==65 && face1==OCT_DIR_RIGHT && face2==OCT_DIR_BACK)
           {
               lookUp=m_uiE2EMapping[elementID*m_uiNumDirections+OCT_DIR_BACK];
               std::cout<<RED<<" element: "<<m_uiAllElements[elementID]<<" RB :"<<m_uiAllElements[lookUp]<<NRM<<std::endl;
           }*/


        /* if(face1==OCT_DIR_LEFT  && face2==OCT_DIR_BACK)
         {
             if(lookUpF1!=LOOK_UP_TABLE_DEFAULT)
                 lookUp=m_uiE2EMapping[lookUpF1*m_uiNumDirections+OCT_DIR_BACK];

         }else if(face1==OCT_DIR_BACK && face2==OCT_DIR_LEFT)
         {
             if(lookUpF1!=LOOK_UP_TABLE_DEFAULT)
                 lookUp=m_uiE2EMapping[lookUpF1*m_uiNumDirections+OCT_DIR_LEFT];
         }else if(face1==OCT_DIR_LEFT && face2==OCT_DIR_FRONT)
         {
             if(lookUpF1!=LOOK_UP_TABLE_DEFAULT && myLev<=m_uiAllElements[lookUpF1].getLevel())
             {
                 lookUp=m_uiE2EMapping[lookUpF1*m_uiNumDirections+OCT_DIR_FRONT];
             }else
             {
                 assert(abs(myLev-m_uiAllElements[lookUpF1].getLevel())==1);


             }
         }*/




        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
        {
            if(myLev>=m_uiAllElements[lookUp].getLevel())
                lookUp=m_uiE2EMapping[lookUp*m_uiNumDirections+face2];
            else
                lookUp=LOOK_UP_TABLE_DEFAULT;

        }else
        {
            lookUp=m_uiE2EMapping[elementID*m_uiNumDirections+face2];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                if(myLev>=m_uiAllElements[lookUp].getLevel())
                    lookUp=m_uiE2EMapping[lookUp*m_uiNumDirections+face1];
                else
                    lookUp=LOOK_UP_TABLE_DEFAULT;
            }

        }



    }



    inline unsigned int Mesh::getDIROfANode(unsigned int ii_x,unsigned int jj_y,unsigned int kk_z)
    {
        // stateX (0-> internal 1-> corner)
        bool stateX[3];
        bool stateY[3];
        bool stateZ[3];


        stateX[0]= ((ii_x>0) && (ii_x<m_uiElementOrder));
        stateX[1]= (ii_x==0);
        stateX[2]= (ii_x==m_uiElementOrder);



        stateY[0]= ((jj_y>0) && (jj_y<m_uiElementOrder));
        stateY[1]= (jj_y==0);
        stateY[2]= (jj_y==m_uiElementOrder);

        stateZ[0]= ((kk_z>0) && (kk_z<m_uiElementOrder));
        stateZ[1]= (kk_z==0);
        stateZ[2]= (kk_z==m_uiElementOrder);

        if(stateX[0] && stateY[0] && stateZ[0])
            return OCT_DIR_INTERNAL;

        if(stateX[1] && stateY[1] && stateZ[1]) // morton index 0
            return OCT_DIR_LEFT_DOWN_BACK;

        if(stateX[2] && stateY[1] && stateZ[1]) // morton index 1
            return OCT_DIR_RIGHT_DOWN_BACK;

        if(stateX[1] && stateY[2] && stateZ[1]) // morton index 2
            return OCT_DIR_LEFT_UP_BACK;

        if(stateX[2] && stateY[2] && stateZ[1]) // morton index 3
            return OCT_DIR_RIGHT_UP_BACK;

        if(stateX[1] && stateY[1] && stateZ[2]) // morton index 4
            return OCT_DIR_LEFT_DOWN_FRONT;

        if(stateX[2] && stateY[1] && stateZ[2]) // morton index 5
            return OCT_DIR_RIGHT_DOWN_FRONT;

        if(stateX[1] && stateY[2] && stateZ[2]) // morton index 6
            return OCT_DIR_LEFT_UP_FRONT;

        if(stateX[2] && stateY[2] && stateZ[2]) // morton index 7
            return OCT_DIR_RIGHT_UP_FRONT;

        // internal faces.

        if(stateZ[0] && stateY[0] && stateX[1])
            return OCT_DIR_LEFT;

        if(stateZ[0] && stateY[0] && stateX[2])
            return OCT_DIR_RIGHT;

        if(stateZ[0] && stateX[0] && stateY[1])
            return OCT_DIR_DOWN;

        if(stateZ[0] && stateX[0] && stateY[2])
            return OCT_DIR_UP;

        if(stateX[0] && stateY[0] && stateZ[1])
            return OCT_DIR_BACK;

        if(stateX[0] && stateY[0] && stateZ[2])
            return OCT_DIR_FRONT;

        // internal corners.

        if(stateZ[0] && stateX[1] && stateY[1])
            return OCT_DIR_LEFT_DOWN;

        if(stateZ[0] && stateX[1] && stateY[2])
            return OCT_DIR_LEFT_UP;

        if(stateY[0] && stateX[1] && stateZ[1])
            return OCT_DIR_LEFT_BACK;

        if(stateY[0] && stateX[1] && stateZ[2])
            return OCT_DIR_LEFT_FRONT;


        if(stateZ[0] && stateX[2] && stateY[1])
            return OCT_DIR_RIGHT_DOWN;

        if(stateZ[0] && stateX[2] && stateY[2])
            return OCT_DIR_RIGHT_UP;

        if(stateY[0] && stateX[2] && stateZ[1])
            return OCT_DIR_RIGHT_BACK;

        if(stateY[0] && stateX[2] && stateZ[2])
            return OCT_DIR_RIGHT_FRONT;



        if(stateX[0] && stateZ[1] && stateY[1])
            return OCT_DIR_DOWN_BACK;

        if(stateX[0] && stateZ[2] && stateY[1])
            return OCT_DIR_DOWN_FRONT;


        if(stateX[0] && stateZ[1] && stateY[2])
            return OCT_DIR_UP_BACK;

        if(stateX[0] && stateZ[2] && stateY[2])
            return OCT_DIR_UP_FRONT;

        assert(false);
        return 0;

    }


}// end of namespace ot.

