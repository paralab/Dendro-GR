/**
 * @file ts.h
 * @author Milinda Fernando. (milinda@cs.utah.edu)
 * @brief Main time stepper information. 
 * @version 0.1
 * @date 2019-12-21
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include<iostream>
#include "dendro.h"
#include "parUtils.h"
#include "dvec.h"
namespace ts
{   

    


    /**
     * @brief default available time integrator methods. 
     */
    enum ETSType {RK3=0, RK4, RK5};
    typedef ot::DVector<DendroScalar, unsigned int> DVec;

    /**@brief data type to store the time stepper level. */
    typedef unsigned int TS_LEV;
    #define ETS_STAGE_DEFAULT 0

    /**@brief : block time stamp*/
    struct BLK_STMP
    {
        TS_LEV _time  = LOOK_UP_TABLE_DEFAULT;
        TS_LEV _stage = ETS_STAGE_DEFAULT ;
        TS_LEV _lev   = LOOK_UP_TABLE_DEFAULT;

        inline bool operator==(const BLK_STMP& other) const {
            return (_time == other._time && _stage==other._stage && _lev ==other._lev);}
        
        inline bool operator!=(const BLK_STMP& other) const {
            return (_time != other._time || _stage!=other._stage || _lev !=other._lev);}
        
    };

    // std::ostream& operator<<(std::ostream& os, BLK_STMP const& other) {
    //         return (os << other._time << " " << other._stage << " " << other._lev);
    // } //end fn

    /**@brief time stepper infomation. */
    struct TSInfo
    {
        /**@brief: time begin*/
        DendroScalar _m_uiTb;

        /**@brief: time end*/
        DendroScalar _m_uiTe;

        /**@brief: current step */
        unsigned int _m_uiStep;
        
        /**@brief: current time*/
        DendroScalar _m_uiT;

        /**@brief: current step size*/
        DendroScalar _m_uiTh;
    };


    /**@brief: Time step flag*/
    enum TSFlag { HALF_STEP = 0, FULL_STEP, NOT_SPECIFIED };

    template<typename T>
    struct BlockTSInfo
    {   
        /**@brief: local block id*/
        T blkID;
        
        /**@brief: time step size flag. */
        TSFlag tsFlg;

        BlockTSInfo(T id, TSFlag ts)
        {
            blkID= id;
            tsFlg = ts;
            
        }

        // void operator==(const BlockTSInfo<T>& other)
        // {
        //     blkID = other.blkID;
        //     tsFlg = other.tsFlg;
        // }

    };

    template<typename T>
    struct SyncBLK
    {
        std::vector<BlockTSInfo<T> > blks;
        /**@brief: lower bound for the level (ghost sync)*/
        unsigned int ll=m_uiMaxDepth;

        

        void append_from_prev(const SyncBLK<T>& lprev)
        {
            for(unsigned int i=0; i< lprev.blks.size(); i++)
                blks.push_back(BlockTSInfo<T>(lprev.blks[i].blkID,TSFlag::FULL_STEP));

            if(ll > lprev.ll)
                ll = lprev.ll;
        }

        // void operator==(const SyncBLK<T>& other)
        // {
        //     blks = other.blks;
        //     ll = other.ll;
        // }

        inline const std::vector<BlockTSInfo<T>>& get_blks() const { return blks;} 

    };


    template<typename T>
    struct SyncBLKWS
    {
        std::vector<T> zVec0;
        std::vector<T> zVec1;

    };



}// end of name space ts. 

namespace par {

    //Forward Declaration
    template <typename T>
    class Mpi_datatype;


    template<>
    class Mpi_datatype<ts::BLK_STMP>
    {
        public :
            static MPI_Datatype value()
            {
                static bool         first = true;
                static MPI_Datatype datatype;

                if (first)
                {
                    first = false;
                    MPI_Type_contiguous(sizeof(ts::BLK_STMP), MPI_BYTE, &datatype);
                    MPI_Type_commit(&datatype);
                }

                return datatype;
            }
    };

}


