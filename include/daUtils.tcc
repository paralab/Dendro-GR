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
        T lagrangeInterpElementToCoord(const ot::Mesh * mesh,const T* in, double * domain_coord, const Point& pt_min, const Point& pt_max , unsigned int elementID,unsigned int interpOrder)
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

            Point d_pt(domain_coord[0],domain_coord[1], domain_coord[2]);

            //!!! Note: This can only happens due to floating point rounding off issues, -------

            if(fabs(domain_coord[0]-pt_min.x()) < 1e-6)
                domain_coord[0] = pt_min.x();

            if(fabs(domain_coord[1]-pt_min.y()) < 1e-6)
                domain_coord[1] = pt_min.y();
                
            if(fabs(domain_coord[2]-pt_min.z()) < 1e-6)
                domain_coord[2] = pt_min.z();

            if(fabs(domain_coord[0]-pt_max.x()) < 1e-6)
                domain_coord[0] = pt_max.x();

            if(fabs(domain_coord[1]-pt_max.y()) < 1e-6)
                domain_coord[1] = pt_max.y();
                
            if(fabs(domain_coord[2]-pt_max.z()) < 1e-6)
                domain_coord[2]= pt_min.z();

            // !!! -------------------------------------------------------------------------------


            if( domain_coord[0] < pt_min.x() || domain_coord[0] > pt_max.x() ) {
                std::cout<<"Error in "<<__func__<<" x coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }
            
            if(domain_coord[1] < pt_min.y() || domain_coord[1] > pt_max.y() ) 
            {
                std::cout<<"Error in "<<__func__<<" y coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }
            
            if(domain_coord[2] < pt_min.z() || domain_coord[2] > pt_max.z() ) 
            {
                std::cout<<"Error in "<<__func__<<" z coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }
            

            
            

            std::vector<T> nodalValues;
            nodalValues.resize(nPe);

            mesh->getElementNodalValues(in,&(*(nodalValues.begin())),elementID);


            const double x=domain_coord[0];
            const double y=domain_coord[1];
            const double z=domain_coord[2];

            const unsigned int shiftSzX = eleOrder / interpOrder;
            const unsigned int shiftSzY = eleOrder / interpOrder;
            const unsigned int shiftSzZ = eleOrder / interpOrder;

            const double eleX=pt_min.x();
            const double eleY=pt_min.y();
            const double eleZ=pt_min.z();

            // step sizes for the interpolation order.
            const double hx = (pt_max.x() - pt_min.x() ) / (interpOrder);
            const double hy = (pt_max.y() - pt_min.y() ) / (interpOrder);
            const double hz = (pt_max.z() - pt_min.z() ) / (interpOrder);

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


        template <typename T>
        T linear_lagrange(const ot::Mesh * mesh,const T* in, double * domain_coord, const Point& pt_min, const Point&  pt_max , unsigned int elementID)
        {
            
            if(!mesh->isActive()) return (T)0.0;

            const unsigned int localBegin=mesh->getElementLocalBegin();
            const unsigned int localEnd=mesh->getElementLocalEnd();

            const unsigned int eleOrder=mesh->getElementOrder();
            const unsigned int nPe=mesh->getNumNodesPerElement();

            const unsigned int rankActive=mesh->getMPIRank();
            const unsigned int npesActive=mesh->getMPICommSize();
            const ot::TreeNode * allElements=&(*(mesh->getAllElements().begin()));


            Point d_pt(domain_coord[0],domain_coord[1], domain_coord[2]);

            //!!! Note: This can only happens due to floating point rounding off issues, -------

            if(fabs(domain_coord[0]-pt_min.x()) < 1e-6)
                domain_coord[0] = pt_min.x();

            if(fabs(domain_coord[1]-pt_min.y()) < 1e-6)
                domain_coord[1] = pt_min.y();
                
            if(fabs(domain_coord[2]-pt_min.z()) < 1e-6)
                domain_coord[2] = pt_min.z();

            if(fabs(domain_coord[0]-pt_max.x()) < 1e-6)
                domain_coord[0] = pt_max.x();

            if(fabs(domain_coord[1]-pt_max.y()) < 1e-6)
                domain_coord[1] = pt_max.y();
                
            if(fabs(domain_coord[2]-pt_max.z()) < 1e-6)
                domain_coord[2]= pt_min.z();

            // !!! -------------------------------------------------------------------------------

            

            if( domain_coord[0] < pt_min.x() || domain_coord[0] > pt_max.x() ) {
                std::cout<<"Error in "<<__func__<<" x coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }
            
            if(domain_coord[1] < pt_min.y() || domain_coord[1] > pt_max.y() ) 
            {
                std::cout<<"Error in "<<__func__<<" y coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }
            
            if(domain_coord[2] < pt_min.z() || domain_coord[2] > pt_max.z() ) 
            {
                std::cout<<"Error in "<<__func__<<" z coord is not lies in the element for Lagrange interpolation with the specified element "<<allElements[elementID]<<std::endl;
                std::cout<<"d_point : "<<d_pt<< " pt_min : "<<pt_min<<" pt_max: "<<pt_max<<std::endl;
                exit(0);
            }

            std::vector<T> nodalValues;
            nodalValues.resize(nPe);

            mesh->getElementNodalValues(in,&(*(nodalValues.begin())),elementID);


            const double x=domain_coord[0];
            const double y=domain_coord[1];
            const double z=domain_coord[2];

            const double eleX=pt_min.x();
            const double eleY=pt_min.y();
            const double eleZ=pt_min.z();

            // step sizes for the interpolation order.
            const double hx = (pt_max.x() - pt_min.x() ) / (eleOrder);
            const double hy = (pt_max.y() - pt_min.y() ) / (eleOrder);
            const double hz = (pt_max.z() - pt_min.z() ) / (eleOrder);

            const unsigned int iix = (unsigned int) std::floor((x - eleX)/hx);
            const unsigned int jjy = (unsigned int) std::floor((y - eleY)/hy);
            const unsigned int kkz = (unsigned int) std::floor((z - eleZ)/hz);

            unsigned int k_index[2] ={ kkz, kkz + 1};
            unsigned int j_index[2] ={ jjy, jjy + 1};
            unsigned int i_index[2] ={ iix, iix + 1};


            assert ( iix < (eleOrder+1) && iix>=0 );
            assert ( jjy < (eleOrder+1) && jjy>=0 );
            assert ( kkz < (eleOrder+1) && kkz>=0 );

            std::vector<T> interp_yz; // x=coord[0] plane (yz) plane interpolated values.
            std::vector<T> interp_z;  // (line x=coord[0], and y=cood[1])

            const unsigned int interp_order = 1;

            interp_yz.resize((interp_order+1)*(interp_order+1),(T)0.0);
            interp_z.resize((interp_order+1),T(0.0));

            /*
                Point simplex_coord[NUM_CHILDREN];

                simplex_coord [ 0 ] = Point( eleX +      iix * hx ,   eleY +  jjy * hy      ,   eleZ + kkz * hz);
                simplex_coord [ 1 ] = Point( eleX +  (iix+1) * hx ,   eleY +  jjy * hy      ,   eleZ + kkz * hz);
                simplex_coord [ 2 ] = Point( eleX +      iix * hx ,   eleY + (jjy + 1) * hy ,   eleZ + kkz * hz);
                simplex_coord [ 3 ] = Point( eleX +  (iix+1) * hx ,   eleY + (jjy + 1) * hy ,   eleZ + kkz * hz);

                simplex_coord [ 4 ] = Point( eleX +      iix * hx ,   eleY +  jjy * hy      ,   eleZ + (kkz+1) * hz);
                simplex_coord [ 5 ] = Point( eleX +  (iix+1) * hx ,   eleY +  jjy * hy      ,   eleZ + (kkz+1) * hz);
                simplex_coord [ 6 ] = Point( eleX +      iix * hx ,   eleY + (jjy + 1) * hy ,   eleZ + (kkz+1) * hz);
                simplex_coord [ 7 ] = Point( eleX +  (iix+1) * hx ,   eleY + (jjy + 1) * hy ,   eleZ + (kkz+1) * hz);
            */


            double Lpoly[2]; // Lagrange polynomials evaluated in in the x direction.
            double sc[2];

            sc[0] = eleX +      (iix) * hx;
            sc[1] = eleX +  (iix + 1) * hx;

            assert(  x >= sc[0] && x <= sc[1] );

            for(unsigned int li = 0 ; li < (interp_order + 1); li+=1)
            {
                Lpoly[li] = 1.0 ;
                for(unsigned int lj = 0 ; lj < (interp_order + 1); lj+=1)
                {
                    if(lj==li) continue;
                    Lpoly[li] *= ((x - sc[lj])) / (sc[li] - sc[lj]);
                }
            }
            

            // -- filling the yz plane.
            for(unsigned int k = 0 ; k < (interp_order + 1); k++)
            for(unsigned int j = 0 ; j < (interp_order + 1); j++)
            {
                interp_yz[k*(interp_order+1)+j]=(T)0.0;
                for(unsigned int li=0; li < (interp_order+1);  li+=1)
                    interp_yz[k*(interp_order+1)+j] +=  nodalValues[k_index[k]*(eleOrder+1)*(eleOrder+1)+j_index[j]*(eleOrder+1) + i_index[li] ] *Lpoly[li];
            }

            sc[0] = eleY +      (jjy) * hy;
            sc[1] = eleY +  (jjy + 1) * hy;

            assert(  y >= sc[0] && y <= sc[1] );


            // -- Lagrange polynomials evaluated in the y direction.
            for(unsigned int li=0;  li<(interp_order+1); li+=1)
            {
                Lpoly[li]=1.0;
                for(unsigned int lj=0; lj<(interp_order+1); lj+=1)
                {
                    if(lj==li) continue;
                    Lpoly[li] *= ((y - sc[lj])) / (sc[li] - sc[lj]);
                }
            }


            for(unsigned int k=0; k<(interp_order+1); k++) 
            {
                interp_z[k]=(T)0.0;
                for(unsigned int li=0; li<interp_order+1; li+=1)
                { 
                    // -- lagrange basis polynomial
                    interp_z[k] += interp_yz[ k * (interp_order+1) + li] * Lpoly[li];
                }
            }

            sc[0] = eleZ +      (kkz) * hz;
            sc[1] = eleZ +  (kkz + 1) * hz;

            assert(  z >= sc[0] && z <= sc[1] );


            // -- Lagrange polynomials evaluated in the z direction.
            for(unsigned int li=0; li < (interp_order+1) ; li+=1)
            {
                Lpoly[li] = 1.0;
                for(unsigned int lj=0; lj<(interp_order+1); lj+=1)
                {
                    if(lj==li) continue;
                    Lpoly[li] *= ((z-sc[lj])) / (sc[li] - sc[lj]);
                }
            }

            T val=(T)0.0;
            for(unsigned int li = 0; li < interp_order + 1; li+=1) 
                val+=interp_z[li]*Lpoly[li];
            
            return val;


        }



        template<typename T, typename CoordT>
        void interpolateToCoords(const ot::Mesh * mesh, const T* in, const CoordT* domain_coords, unsigned int length, const Point* const grid_limit, const Point* const domain_limit, T* out,std::vector<unsigned int >& validIndices)
        {

            unsigned int rankGlobal=mesh->getMPIRankGlobal();
            unsigned int npesGlobal=mesh->getMPICommSizeGlobal();

            if(mesh->isActive())
            {
                unsigned int rankActive=mesh->getMPIRank();
                unsigned int npesActive=mesh->getMPICommSize();

                const double gridRangeX = grid_limit[1].x() - grid_limit[0].x(); // x range of the grid.
                const double gridRangeY = grid_limit[1].y() - grid_limit[0].y(); // y range of the grid.
                const double gridRangeZ = grid_limit[1].z() - grid_limit[0].z(); // z range of the grid.

                const double domainRangeX = domain_limit[1].x() - domain_limit[0].x(); // x range of the domain.
                const double domainRangeY = domain_limit[1].y() - domain_limit[0].y(); // y range of the domain. 
                const double domainRangeZ = domain_limit[1].z() - domain_limit[0].z(); // z range of the domain. 

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
                unsigned int octree_coords[m_uiDim];
                
                for(unsigned int i=0;i<numPts;i++)
                {   
                    octree_coords[ 0 ] = (unsigned int) (grid_limit[0].x() + (domain_coords[m_uiDim*i + 0] - domain_limit[0].x())*(gridRangeX/domainRangeX));
                    octree_coords[ 1 ] = (unsigned int) (grid_limit[0].y() + (domain_coords[m_uiDim*i + 1] - domain_limit[0].y())*(gridRangeY/domainRangeY));
                    octree_coords[ 2 ] = (unsigned int) (grid_limit[0].z() + (domain_coords[m_uiDim*i + 2] - domain_limit[0].z())*(gridRangeZ/domainRangeZ));

                    coordOcts_skey.push_back(ot::SearchKey(octree_coords[0], octree_coords[1], octree_coords[2], m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
                    coordOcts_skey.back().addOwner(i);
                }

                // 2. get the element local splitters.
                const std::vector<ot::TreeNode>& splitterElements=mesh->getSplitterElements();
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
                Point pt_min;
                Point pt_max;

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

                        coord[0] = domain_limit[0].x()  +  ((meshOctree[searchResult].minX() - grid_limit[0].x())*(domainRangeX/gridRangeX)) ;
                        coord[1] = domain_limit[0].y()  +  ((meshOctree[searchResult].minY() - grid_limit[0].y())*(domainRangeY/gridRangeY)) ;
                        coord[2] = domain_limit[0].z()  +  ((meshOctree[searchResult].minZ() - grid_limit[0].z())*(domainRangeZ/gridRangeZ)) ;

                        pt_min = Point(coord[0], coord[1], coord[2]);

                        coord[0] = domain_limit[0].x()  +  ((meshOctree[searchResult].maxX() - grid_limit[0].x())*(domainRangeX/gridRangeX)) ;
                        coord[1] = domain_limit[0].y()  +  ((meshOctree[searchResult].maxY() - grid_limit[0].y())*(domainRangeY/gridRangeY)) ;
                        coord[2] = domain_limit[0].z()  +  ((meshOctree[searchResult].maxZ() - grid_limit[0].z())*(domainRangeZ/gridRangeZ)) ;

                        pt_max = Point(coord[0], coord[1], coord[2]);
                        
                        // perform the linterpolation
                        for(unsigned int w=0;w<ownerList->size();w++)
                        {
                            coord[0]=domain_coords[m_uiDim * (*ownerList)[w]     ];
                            coord[1]=domain_coords[m_uiDim * (*ownerList)[w] + 1 ];
                            coord[2]=domain_coords[m_uiDim * (*ownerList)[w] + 2 ];
                            out[(*ownerList)[w]]=lagrangeInterpElementToCoord(mesh,in,coord, pt_min, pt_max, searchResult,mesh->getElementOrder());
                            //out[(*ownerList)[w]] = linear_lagrange(mesh, in, coord, pt_min, pt_max, searchResult);
                            validIndices.push_back((*ownerList)[w]);
                        }

                    }
                
                }

                std::sort(validIndices.begin(),validIndices.end());
                meshOctree.clear();


            }


        }

        template<typename T, typename CoordT>
        void interpolateToCoordsAndGather(const ot::Mesh * mesh, const T* in, const CoordT* domain_coords, unsigned int length, const Point* const grid_limit, const Point* const domain_limit, T* out,unsigned int root,unsigned int dof)
        {

            if(mesh->isActive())
            {
                MPI_Comm comm_active = mesh->getMPICommunicator();
                unsigned int rank_active = mesh->getMPIRank();
                unsigned int npes_active = mesh->getMPICommSize();

                std::vector<unsigned int > v_index;
                v_index.reserve(length);

                unsigned int in_sz  = mesh->getDegOfFreedom();
                unsigned int out_sz = length; 
                
                interpolateToCoords(mesh,in + 0*in_sz,domain_coords,length,grid_limit,domain_limit,out + 0 * out_sz,v_index);
                int s_count = v_index.size() * dof;
                std::vector<int> r_count;
                std::vector<int> r_offset;
                
                if(rank_active == root)
                {
                    r_count.resize(npes_active);
                    r_offset.resize(npes_active);
                }
                    
                par::Mpi_Gather(&s_count,r_count.data(),1,root,comm_active);

                if(rank_active == root)
                {
                    r_offset[0]=0;
                    omp_par::scan(r_count.data(),r_offset.data(),r_count.size());
                }
                    

                for(unsigned int v_id =1; v_id < dof; v_id++)
                {
                    v_index.clear();
                    interpolateToCoords(mesh,in + v_id*in_sz,domain_coords,length,grid_limit,domain_limit,out + v_id * out_sz,v_index);
                }

                std::vector<T> send_buff;
                send_buff.resize(s_count);

                for(unsigned int w=0; w< v_index.size(); w++)
                    for(unsigned int v_id =0; v_id < dof; v_id++)
                        send_buff[w * dof + v_id] = out[v_id * out_sz + v_index[w]];
                

                par::Mpi_Gatherv(send_buff.data(),s_count,out,r_count.data(),r_offset.data(),root,comm_active);



            }
            return;

        }

    } // end of namespace da

}// end of namespace ot