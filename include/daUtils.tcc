//
// Created by milinda on 4/27/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah.
 * @brief contains utilitity functions of mesh based distributed array(da)
 * */
#include "mesh.h"
#include "TreeNode.h"
#include "sfcSort.h"
#include "sfcSearch.h"
#include <iostream>
namespace ot
{

    namespace da
    {


        template <typename T>
        T lagrangeInterpElementToCoord(const ot::Mesh * mesh,const T* in,double * coord,unsigned int elementID,unsigned int interpOrder)
        {

            if(!mesh->isActive()) return (T)0.0;

            const unsigned int localBegin=mesh->getElementLocalBegin();
            const unsigned int localEnd=mesh->getElementLocalEnd();

            const unsigned int eleOrder=mesh->getElementOrder();
            const unsigned int nPe=mesh->getNumNodesPerElement();

            const unsigned int rankActive=mesh->getMPIRank();
            const unsigned int npesActive=mesh->getMPICommSize();
            const ot::TreeNode * allElements=&(*(mesh->getAllElements().begin()));

            if(interpOrder>eleOrder) {std::cout<<" Error in "<<__func__<<" interpOrder "<<interpOrder<<" > eleOrder "<<eleOrder<<std::endl; return (T)0.0;}
            if(eleOrder%interpOrder!=0) {std::cout<<" Error in "<<__func__<<" eleOrder "<<eleOrder<<" not perfectly divisible by interpOrder "<<interpOrder<<std::endl; return (T)0.0;}

            if(coord[0]<allElements[elementID].minX() || coord[0]>allElements[elementID].maxX() ) {std::cout<<"Error in "<<__func__<<" x coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;}
            if(coord[1]<allElements[elementID].minY() || coord[1]>allElements[elementID].maxY() ) {std::cout<<"Error in "<<__func__<<" y coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;}
            if(coord[2]<allElements[elementID].minZ() || coord[2]>allElements[elementID].maxZ() ) {std::cout<<"Error in "<<__func__<<" z coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;}


            std::vector<T> nodalValues;
            nodalValues.resize(nPe);

            mesh->getElementNodalValues(in,&(*(nodalValues.begin())),elementID);


            const double x=coord[0];
            const double y=coord[1];
            const double z=coord[2];

            const unsigned int shiftSzX=eleOrder/interpOrder;
            const unsigned int shiftSzY=eleOrder/interpOrder;
            const unsigned int shiftSzZ=eleOrder/interpOrder;

            const double eleX=(double)allElements[elementID].getX();
            const double eleY=(double)allElements[elementID].getY();
            const double eleZ=(double)allElements[elementID].getZ();

            // step sizes for the interpolation order.
            const double hx=(double)(allElements[elementID].maxX()-allElements[elementID].minX())/(interpOrder);
            const double hy=(double)(allElements[elementID].maxY()-allElements[elementID].minY())/(interpOrder);
            const double hz=(double)(allElements[elementID].maxZ()-allElements[elementID].minZ())/(interpOrder);

            std::vector<T> interp_yz; // x=coord[0] plane (yz) plane interpolated values.
            std::vector<T> interp_z;  // (line x=coord[0], and y=cood[1])
            interp_yz.resize((eleOrder+1)*(eleOrder+1),(T)0.0);
            interp_z.resize(eleOrder+1,T(0.0));

            // computation of the lagrange polynomials in the x direction.
            double Lpoly[interpOrder+1]; // Lagrange polynomials evaluated in in the x direction.
            for(unsigned int li=0;li<(eleOrder+1);li+=shiftSzX)
            {
                Lpoly[li/shiftSzX]=1.0;
                for(unsigned int lj=0;lj<(eleOrder+1);lj+=shiftSzX)
                {
                    if(lj==li) continue;
                    Lpoly[li/shiftSzX]*=((x-(eleX+hx*lj)))/((eleX+hx*li)-(eleX+hx*lj));
                }
            }


            // filling the yz plane.
            for(unsigned int k=0;k<(eleOrder+1);k++)
                for(unsigned int j=0;j<(eleOrder+1);j++)
                {
                    interp_yz[k*(eleOrder+1)+j]=(T)0.0;
                    for(unsigned int li=0;li<(eleOrder+1);li+=shiftSzX)
                        interp_yz[k*(eleOrder+1)+j]+= nodalValues[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+li]*Lpoly[li/shiftSzX];
                }


            // Lagrange polynomials evaluated in the y direction.
            for(unsigned int li=0;li<(eleOrder+1);li+=shiftSzY)
            {
                Lpoly[li/shiftSzY]=1.0;
                for(unsigned int lj=0;lj<(eleOrder+1);lj+=shiftSzY)
                {
                    if(lj==li) continue;
                    Lpoly[li/shiftSzY]*=((y-(eleY+hy*lj)))/((eleY+hy*li)-(eleY+hy*lj));
                }
            }


            for(unsigned int k=0;k<(eleOrder+1);k++) // (line x=coord[0], and y=cood[1])
            {
                interp_z[k]=(T)0.0;
                for(unsigned int li=0;li<eleOrder+1;li+=shiftSzY)
                { // lagrange basis polynomial
                    interp_z[k]+=interp_yz[k*(eleOrder+1)+li]*Lpoly[li/shiftSzY];
                }
            }


            // Lagrange polynomials evaluated in the z direction.
            for(unsigned int li=0;li<(eleOrder+1);li+=shiftSzZ)
            {
                Lpoly[li/shiftSzZ]=1.0;
                for(unsigned int lj=0;lj<(eleOrder+1);lj+=shiftSzZ)
                {
                    if(lj==li) continue;
                    Lpoly[li/shiftSzY]*=((z-(eleZ+hz*lj)))/((eleZ+hz*li)-(eleZ+hz*lj));
                }
            }

            T val=(T)0.0;
            for(unsigned int li=0;li<eleOrder+1;li+=shiftSzX) // interpolating to the point (coord[0],coord[1],coord[2])
            { // lagrange basis polynomial
               val+=interp_z[li]*Lpoly[li/shiftSzZ];
            }

            return val;

        }



        template<typename T, typename CoordT>
        void interpolateToCoords(const ot::Mesh * mesh, const T* in, const CoordT* coords,unsigned int length,T* out,std::vector<unsigned int >& validIndices)
        {

            unsigned int rankGlobal=mesh->getMPIRankGlobal();
            unsigned int npesGlobal=mesh->getMPICommSizeGlobal();

            if(mesh->isActive())
            {
                unsigned int rankActive=mesh->getMPIRank();
                unsigned int npesActive=mesh->getMPICommSize();

                const ot::TreeNode* allElements=&(*(mesh->getAllElements().begin()));

                std::vector<ot::TreeNode> meshOctree;
                meshOctree.resize(mesh->getAllElements().size());

                for(unsigned int ele=0;ele<meshOctree.size();ele++)
                    meshOctree[ele]=allElements[ele];

                const unsigned int localElementBegin=mesh->getElementLocalBegin();
                const unsigned int localElementEnd=mesh->getElementLocalEnd();

                const unsigned int eleOrder=mesh->getElementOrder();
                const unsigned int nPe=mesh->getNumNodesPerElement();

                std::vector<T> nodalValues;
                nodalValues.resize(nPe);

                if(length % m_uiDim != 0 )
                {
                    std::cout<<"Invalid use of "<<__func__<<" length of coords should be perfectly divisible by m_uiDim: "<<m_uiDim<<std::endl;
                    return ;
                }

                // 1. convert coords to octants
                const unsigned int numPts=length/m_uiDim;
                std::vector<ot::SearchKey> coordOcts_skey;
                for(unsigned int i=0;i<numPts;i++)
                {
                    coordOcts_skey.push_back(ot::SearchKey((unsigned int)(coords[m_uiDim*i]),(unsigned int)(coords[m_uiDim*i+1]),(unsigned int)(coords[m_uiDim*i+2]),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
                    coordOcts_skey.back().addOwner(i);
                }

                // 2. get the element local splitters.
                const std::vector<ot::TreeNode> splitterElements=mesh->getSplitterElements();
                std::vector<ot::Key> splitterKeys;  // contains the splitter key begin.
                for(unsigned int i=0;i<npesActive;i++)
                {
                    coordOcts_skey.push_back(ot::SearchKey(splitterElements[2*i]));
                    splitterKeys.push_back(ot::Key(splitterElements[2*i]));
                }

                // 3. merge the coords keys.
                std::vector<ot::Key> coordOcts_key;
                mergeKeys(coordOcts_skey,coordOcts_key); // note that this will result sorted keys.
                coordOcts_skey.clear();

                // 4. Search the splitter keys in the coordOCt_keys to find which coords resides in the local proc.
                SFC::seqSearch::SFC_treeSearch(&(*(splitterKeys.begin())),&(*(coordOcts_key.begin())),0,splitterKeys.size(),0,coordOcts_key.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

                assert(splitterKeys[rankActive].getFlag() & OCT_FOUND);
                if(!(splitterKeys[rankActive].getFlag() & OCT_FOUND))
                {
                    std::cout<<"Error["<<__func__<<"]: splitter keys not found "<<std::endl;
                    return;
                }

                const unsigned int myBegin=splitterKeys[rankActive].getSearchResult();
                const unsigned int myEnd = (rankActive+1<npesActive) ?  splitterKeys[rankActive+1].getSearchResult() : coordOcts_key.size();

                validIndices.clear();
                unsigned int outSz=myEnd-myBegin;

                // 5. Search the coordOcts_key in the allElements to find out the corresponding keys.
                SFC::seqSearch::SFC_treeSearch(&(*(coordOcts_key.begin())),&(*(meshOctree.begin())),myBegin,myEnd,localElementBegin,localElementEnd,m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

                unsigned int searchResult;
                double coord[3];
                std::vector<unsigned int >* ownerList;
                for(unsigned int i=myBegin;i<myEnd;i++)
                {
                    assert(coordOcts_key[i].getFlag() & OCT_FOUND);
                    if(!(coordOcts_key[i].getFlag() & OCT_FOUND))
                    {
                        std::cout<<"rank: "<<rankActive<<" key: "<<coordOcts_key[i]<<"not found "<<std::endl;
                    }
                    if(coordOcts_key[i].getFlag() & OCT_FOUND)
                    {
                        searchResult=coordOcts_key[i].getSearchResult();
                        ownerList=coordOcts_key[i].getOwnerList();
                        // perform the linterpolation
                        for(unsigned int w=0;w<ownerList->size();w++)
                        {
                            coord[0]=coords[3*(*ownerList)[w]];
                            coord[1]=coords[3*(*ownerList)[w]+1];
                            coord[2]=coords[3*(*ownerList)[w]+2];
                            out[(*ownerList)[w]]=lagrangeInterpElementToCoord(mesh,in,coord,searchResult,4);
                            validIndices.push_back((*ownerList)[w]);
                        }

                    }
                }

                std::sort(validIndices.begin(),validIndices.end());

                meshOctree.clear();


            }


        }

    } // end of namespace da

}// end of namespace ot