/**
 * @file waveletRefEl.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Wavelet reference element class. 
 * @version 0.1
 * @date 2020-06-10
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#include "waveletRefEl.h"

namespace wavelet
{
    WaveletEl::WaveletEl(RefElement* refEl)
    {
        m_uiRefEl = refEl;
        const unsigned int eleOrder = m_uiRefEl->getElementOrder();
        assert(eleOrder%2==0); // currently coded for even element order. 
        const unsigned int pw       = eleOrder>>1u;

        m_uiRefElBdy = new RefElement(1,eleOrder/2);

        const unsigned int wb = eleOrder>>1u;
        const unsigned int we = eleOrder + wb+1;

        // for(unsigned int i=wb-1; i < we; i+=2)
        //     m_uiCIndex.push_back(i);

        const unsigned int ib = 0; 
        const unsigned int ie = eleOrder + 2*pw +1;
        
        for(unsigned int i=ib;i < ie; i+=2)
            m_uiPIndex.push_back(i);

        for(unsigned int i=ib+1; i < ie; i+=2)
            m_uiCIndex.push_back(i);


        const double * I0 = m_uiRefEl->getIMTChild0();
        const double * I1 = m_uiRefEl->getIMTChild1();

        const unsigned int num_wc_1d = m_uiCIndex.size();
        // m_uiWInterpOp.resize(num_wc_1d*(eleOrder+1));
        // double * wOp      = m_uiWInterpOp.data();

        
        // unsigned int wcount = 0; 
        // for(unsigned int cnum=0; cnum < 2; cnum++)
        // {
        //     for(unsigned int m=0; m < eleOrder+1 ; m++)
        //     {
        //         const unsigned int cid = (cnum * eleOrder) + m;
        //         if(cid == m_uiCIndex[wcount] )
        //         {
        //             if(cnum==0)
        //              std::memcpy(wOp + wcount*(eleOrder+1) , I0 + m*(eleOrder+1), sizeof(double)*(eleOrder+1));
        //             else
        //              std::memcpy(wOp + wcount*(eleOrder+1) , I1 + m*(eleOrder+1), sizeof(double)*(eleOrder+1));
        //             wcount++;
        //         }
        //     }
        // }

        // std::cout<<"I0:"<<std::endl;
        // printArray_2D(I0,eleOrder+1,eleOrder+1);
        
        // std::cout<<"I1:"<<std::endl;
        // printArray_2D(I1,eleOrder+1,eleOrder+1);

        //std::cout<<"pIndex:"<<std::endl;
        //printArray_1D(m_uiPIndex.data(),m_uiPIndex.size());
        
        //std::cout<<"cIndex:"<<std::endl;
        //printArray_1D(m_uiCIndex.data(),m_uiCIndex.size());


        const unsigned int nPe = (eleOrder+1)*(eleOrder+1)*(eleOrder+1);
        const unsigned int nPe_children = (2*eleOrder+1)*(2*eleOrder+1)*(2*eleOrder+1);

        m_uiVIn.resize(nPe);
        m_uiNVec.resize(nPe);
        m_uiVOut.resize(nPe_children);

    }

    WaveletEl::~WaveletEl()
    {

        delete m_uiRefElBdy;
        m_uiNVec.clear();
        m_uiVIn.clear();
        m_uiVOut.clear();
    }

    void WaveletEl::compute_wavelets_3D(const double* in, const unsigned int * isz ,std::vector<double>& wc, bool isBdy) 
    {
        
        const unsigned int inx = isz[0];
        const unsigned int iny = isz[1];
        const unsigned int inz = isz[2];

        const unsigned int eleOrder = m_uiRefEl->getElementOrder();
        const unsigned int np_1d= eleOrder+1;
        const unsigned int nPe = np_1d*np_1d*np_1d;

        const unsigned int num_wc_1d = m_uiCIndex.size();
        const unsigned int num_wc = num_wc_1d*num_wc_1d*num_wc_1d;

        const double in_min = 1.0;
        const bool compute_relative_wc= true;

        unsigned int wcount=0;
        wc.resize(num_wc);

        // initialize the wavelets to zero. 
        for(unsigned int i=0; i < wc.size(); i++)
            wc[i]=0;


        if(isBdy)
        {
            // element is a boundary element. 
            const unsigned int bdyEOrder = eleOrder>>1u;
            const unsigned int pw = eleOrder>>1u;

            const unsigned int bnp_1d = (bdyEOrder+1);
            const unsigned int bnp_child_1d = 2*bdyEOrder+1;

            for(unsigned int k=pw; k < (pw + eleOrder+1); k+=2)
            for(unsigned int j=pw; j < (pw + eleOrder+1); j+=2)
            for(unsigned int i=pw; i < (pw + eleOrder+1); i+=2)
            {
                const unsigned int pp  = k*(iny*inx) + j*(inx) + i;
                const unsigned int pp1 = ((k-pw)>>1u) * ( bnp_1d*bnp_1d ) +  ((j-pw)>>1u)*bnp_1d + ((i-pw)>>1u);
                m_uiVIn[pp1] = in[pp];
            }


            int bit[3];
            for(unsigned int cnum=0; cnum < NUM_CHILDREN; cnum++)
            {
                m_uiRefElBdy->I3D_Parent2Child(m_uiVIn.data(),m_uiNVec.data(),cnum);
                
                bit[0] = binOp::getBit(cnum,0);
                bit[1] = binOp::getBit(cnum,1);
                bit[2] = binOp::getBit(cnum,2);

                const unsigned int ib = bdyEOrder*bit[0];
                const unsigned int jb = bdyEOrder*bit[1];
                const unsigned int kb = bdyEOrder*bit[2];

                for(unsigned int k=kb; k < (kb + (bdyEOrder+1)); k++)
                for(unsigned int j=jb; j < (jb + (bdyEOrder+1)); j++)
                for(unsigned int i=ib; i < (ib + (bdyEOrder+1)); i++)
                    m_uiVOut[k*bnp_child_1d*bnp_child_1d + j*bnp_child_1d + i] = m_uiNVec[(k-kb)*bnp_1d*bnp_1d + (j-jb)*bnp_1d + (i-ib)];

            }

            
            if(compute_relative_wc)
            {
                wcount=0;
                for(unsigned int k=pw+1; k < (pw + eleOrder+1); k+=2)
                for(unsigned int j=pw+1; j < (pw + eleOrder+1); j+=2)
                for(unsigned int i=pw+1; i < (pw + eleOrder+1); i+=2,wcount++)
                {
                    const unsigned int pp  = k*(iny*inx) + j*(inx) + i;
                    const unsigned int pp1 = (k-pw) * ( bnp_child_1d*bnp_child_1d ) +  (j-pw) * bnp_child_1d + (i-pw);
                    wc[wcount] = fabs(in[pp]- m_uiVOut[pp1])/std::max( in_min, fabs(in[pp]) ); // relative wavelet tol. 
                }

            }else
            {
                wcount=0;
                for(unsigned int k=pw+1; k < (pw + eleOrder+1); k+=2)
                for(unsigned int j=pw+1; j < (pw + eleOrder+1); j+=2)
                for(unsigned int i=pw+1; i < (pw + eleOrder+1); i+=2,wcount++)
                {
                    const unsigned int pp  = k*(iny*inx) + j*(inx) + i;
                    const unsigned int pp1 = (k-pw) * ( bnp_child_1d*bnp_child_1d ) +  (j-pw) * bnp_child_1d + (i-pw);
                    wc[wcount] = fabs(in[pp]- m_uiVOut[pp1]); // absolute wavelet tol. 
                }

            }



            

            return;


        }

        


        for(unsigned int k=0; k < m_uiPIndex.size(); k++)
        for(unsigned int j=0; j < m_uiPIndex.size(); j++)
        for(unsigned int i=0; i < m_uiPIndex.size(); i++)
        {
            const unsigned int pp = m_uiPIndex[k]*(iny*inx) + m_uiPIndex[j]*(inx) + m_uiPIndex[i];
            m_uiVIn[k*(np_1d)*(np_1d) + j*(np_1d) + i] = in[pp]; 
        }

        
        int bit[3];
        for(unsigned int cnum=0; cnum < NUM_CHILDREN; cnum++)
        {
            m_uiRefEl->I3D_Parent2Child(m_uiVIn.data(),m_uiNVec.data(),cnum);
            
            bit[0] = binOp::getBit(cnum,0);
            bit[1] = binOp::getBit(cnum,1);
            bit[2] = binOp::getBit(cnum,2);

            const unsigned int ib = eleOrder*bit[0];
            const unsigned int jb = eleOrder*bit[1];
            const unsigned int kb = eleOrder*bit[2];

            for(unsigned int k=kb; k < (kb + (eleOrder+1)); k++)
            for(unsigned int j=jb; j < (jb + (eleOrder+1)); j++)
            for(unsigned int i=ib; i < (ib + (eleOrder+1)); i++)
                m_uiVOut[k*iny*inx + j*inx + i] = m_uiNVec[(k-kb)*np_1d*np_1d + (j-jb)*np_1d + (i-ib)];

        }

        // for(unsigned int w=0; w < inx*iny*inz; w++)
        //     std::cout<<"in["<<w<<"]:"<<in[w] <<" \t | \t out["<<w<<"]:"<<m_uiVOut[w]<<std::endl;

        
        // for(unsigned int k=0; k < inz; k++)
        // for(unsigned int j=0; j < iny; j++)
        // for(unsigned int i=0; i < inx; i++)
        // {
        //     unsigned int pp = k*iny*inx + j* inx + i;
        //     printf("in (%d, %d, %d) : %f    out(%d, %d, %d ): %f \n",i,j,k,in[pp],i,j,k,m_uiVOut[pp]);
        // }

        // {
        //     double alpha=1.0;
        //     double beta=0.0;
        //     int incx=1;
        //     int incy=1;
        //     dgemv_((char*)"N",(int*)&num_wc,(int*)&nPe,&alpha,m_uiWInterp_3d.data(),(int*)&num_wc,m_uiVIn.data(),&incx,&beta,m_uiVOut.data(),&incy);

        // }

        
        

        if(compute_relative_wc)
        {
            wcount=0;
            for(unsigned int k=0; k < m_uiCIndex.size(); k++)
            for(unsigned int j=0; j < m_uiCIndex.size(); j++)
            for(unsigned int i=0; i < m_uiCIndex.size(); i++,wcount++)
            {
                const unsigned int pp = m_uiCIndex[k]*(iny*inx) + m_uiCIndex[j]*(inx) + m_uiCIndex[i];
                wc[wcount] = fabs(m_uiVOut[pp] - in[pp])/std::max( in_min, fabs(in[pp]) );

            }

        }else
        {
            wcount=0;
            for(unsigned int k=0; k < m_uiCIndex.size(); k++)
            for(unsigned int j=0; j < m_uiCIndex.size(); j++)
            for(unsigned int i=0; i < m_uiCIndex.size(); i++,wcount++)
            {
                const unsigned int pp = m_uiCIndex[k]*(iny*inx) + m_uiCIndex[j]*(inx) + m_uiCIndex[i];
                wc[wcount] = fabs(m_uiVOut[pp] - in[pp]); // abs wavelets. 
             
            }

        }
        
        
            

    }



}// end of wavelets