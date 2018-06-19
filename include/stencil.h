//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains the implementation of the structure stencil which is generic enough to specify,a given
 * stencil.
*/
//

#ifndef SFCSORTBENCH_STENCIL_H
#define SFCSORTBENCH_STENCIL_H

#include <iostream>
#include <array>
#include "assert.h"

/**
 * @author Milinda Fernando
 * @brief Contains attribute to specify a generic stencil
 * */

enum StencilDirection{STENCIL_DIR_X,STENCIL_DIR_Y,STENCIL_DIR_Z,STENCIL_DIR_XY,STENCIL_DIR_YZ,STENCIL_DIR_XZ,STENCIL_DIR_XYZ};

template  <typename T,unsigned int length, unsigned int offset>
struct Stencil
{

private:
    /**variable to store stencil coefficients*/
    std::array<T,length> m_uiS;
    /**which direction to apply the stencil. */
    StencilDirection m_uiDirection;
    /** offset location of the stencil , which denotes the 0th point coefficient. */
    unsigned int m_uiOffset=offset;


public:
    /** @brief: default constructor */
    Stencil(){};

    /**
     * @brief: creates a stencil with given coefficients.
     * @param[in] coeff: coefficients for the stencil.
     * @param[in] n : number of coefficients.
     */
    Stencil(T* coeff, int n, StencilDirection pDir)
    {
        assert(n==length);
        for(unsigned int k=0;k<length;k++)
            m_uiS[k]=coeff[k];

        m_uiDirection=pDir;
    }

    const T& operator[](unsigned int i) const
    {
        assert(i<length);
        return m_uiS[i];
    }

    /** returns the stencil direction*/
    inline StencilDirection getStencilDirection() const {return m_uiDirection;}
    /** returns the offset of the stencil.*/
    inline unsigned int getOffset() const {return m_uiOffset;}
    /** returns the length of the stencil.*/
    inline unsigned int getStencilLength() const {return m_uiS.size();}


};


#endif //SFCSORTBENCH_STENCIL_H
