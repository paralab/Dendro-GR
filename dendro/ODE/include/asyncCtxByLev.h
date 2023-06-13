/**
 * @file actx.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Asynchronous communication for non uniform time stepper. This class is derived from the asyncExchangeContext.h to support partial ghost sync.
 * @version 0.1
 * @date 2020-01-16
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */

#pragma once

#include "asyncExchangeContex.h"
namespace ot
{
    class AsyncCommByLev : public ot::AsyncExchangeContex
    {
        protected:
            int m_uiLev;
        
        public:
            /**@brief: default constructor*/
            AsyncCommByLev(const void* var, int lev) : ot::AsyncExchangeContex(var)
            {
                m_uiLev = lev;
            }

            /**@brief: overloaded operator for partial ghost sync*/
            bool operator==(const AsyncCommByLev& other)
            {
                return (this->m_uiBuffer == other.m_uiBuffer && m_uiLev == other.m_uiLev);
            }

            /**@brief: returns the level of the ghost exchange. */
            inline unsigned int getLevel() const {return m_uiLev;}

    };

} // end of name space ts