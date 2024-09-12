#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "derivs.h"
#include "fluidUtils.h"
#include "mathUtils.h"
#include "parameters.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

#define deriv_x deriv42_x
#define deriv_y deriv42_y
#define deriv_z deriv42_z

#define deriv_xx deriv42_xx
#define deriv_yy deriv42_yy
#define deriv_zz deriv42_zz

#define adv_deriv_x deriv42adv_x
#define adv_deriv_y deriv42adv_y
#define adv_deriv_z deriv42adv_z

#define ko_deriv_x ko_deriv42_x
#define ko_deriv_y ko_deriv42_y
#define ko_deriv_z ko_deriv42_z

namespace fluid {
/***
 *@brief computes the fluid rhs in unzip format.
 * @param [in] pMesh: current mesh
 * @param [in] time: rk time
 * @param [out] unzipVarsRHS: unzip variables (rhs alloc)
 * @param [in] uZipVars: evolving vars in unzip format.
 * @param [in] uZipVarsPrim: unzip variables for primitive vars.
 * */
void fluid_rhs(const ot::Mesh *pMesh, DendroScalar time,
               DendroScalar **uzipVarsRHS, const DendroScalar **uZipVars,
               const DendroScalar **uZipVarsPrim);

/***
 *@brief  enforce coundary conditions for fluid eqs (evolv vars).
 * @param [in] pMesh: current mesh
 * @param [in/out] uZipVars: evolving vars in unzip format.
 * */
void fluid_bcs_cons(const ot::Mesh *pMesh, DendroScalar **uZipVars);

/***
 *@brief  enforce coundary conditions for fluid eqs (primitive vars).
 * @param [in] pMesh: current mesh
 * @param [in/out] uZipVarsPrim: unzip variables for primitive vars.
 * */
void fluid_bcs_prim(const ot::Mesh *pMesh, DendroScalar **uZipVarsPrim);

/***
 *@brief compute rhs for a single block.
 * @param [in] time: rk time
 * @param [out] unzipVarsRHS: unzip variables (rhs alloc)
 * @param [in] uZipVars: evolving vars in unzip format.
 * @param [in] uZipVarsPrim: unzip variables for primitive vars.
 * @param [in] offset : offset for the unzip representation
 * @param [in] ptmin: point min
 * @param [in] ptmax: point max
 * @param [in] sz: block size
 * @param [in] bflag: boundary flag
 * */
void fluid_block_rhs(DendroScalar time, DendroScalar **uzipVarsRHS,
                     const DendroScalar **uZipVars,
                     const DendroScalar **uZipVarsPrim,
                     const unsigned int &offset, const DendroScalar *ptmin,
                     const DendroScalar *ptmax, const unsigned int *sz,
                     const unsigned int &bflag);

/***
 *@brief enforce the boundary condition for a sinlgle block.
 * @param [in] uZipVars: evolving vars in unzip format.
 * @param [in] offset : offset for the unzip representation
 * @param [in] ptmin: point min
 * @param [in] ptmax: point max
 * @param [in] sz: block size
 * @param [in] bflag: boundary flag
 * */
void fluid_block_bcs_cons(DendroScalar **uZipVars, const unsigned int &offset,
                          const DendroScalar *pmin, const DendroScalar *pmax,
                          const unsigned int *sz, const unsigned int &bflag);

/***
 *@brief enforce the boundary condition for a sinlgle block.([prim vars)
 * @param [in/out] uZipVarsPrim: unzip variables for primitive vars.
 * @param [in] offset : offset for the unzip representation
 * @param [in] ptmin: point min
 * @param [in] ptmax: point max
 * @param [in] sz: block size
 * @param [in] bflag: boundary flag
 * */
void fluid_block_bcs_prim(DendroScalar **uZipVarsPrim,
                          const unsigned int &offset, const DendroScalar *pmin,
                          const DendroScalar *pmax, const unsigned int *sz,
                          const unsigned int &bflag);

void fake_initial_data(DendroScalar x, DendroScalar y, DendroScalar z,
                       DendroScalar *u);

}  // end of namespace fluid

#endif
