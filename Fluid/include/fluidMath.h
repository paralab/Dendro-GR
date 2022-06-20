#ifndef FLUID_MATH_H
#define FLUID_MATH_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <set>
#include "geom.h"

#include "fluid.h"
#include "fluidUtils.h"
#include "mathUtils.h"
#include "parameters.h"
#include <execinfo.h>
#include <cxxabi.h>

// Work array for con_to_prim
#define NFPAR 24

#define FP_SSQ 0
#define FP_BSQ 1
#define FP_BS 2
#define FP_D 3
#define FP_TAU 4
#define FP_GAMMA 5
#define FP_VSQ 6
#define FP_COORD_X 7
#define FP_COORD_Y 8
#define FP_COORD_Z 9
#define FP_CHR 10
#define FP_INDEX_X 11
#define FP_INDEX_Y 12
#define FP_INDEX_Z 13
#define FP_KAPPA 14
#define FP_TRACE 15
#define FP_CAL_PRESS 16
#define FP_Q 17
#define FP_C2PWARN 18
#define FP_VACUUM 19

#define VELOCITY_FACTOR 0.9999

namespace fluid_math
{
/**
 * @brief Error codes and their meaning for the fluid solvers (i.e. cons to prim and vice versa. )
 */
enum Error
{
    FLUID_SOLVER_SUCCESS = 0,           //solver  successful
    FLUID_W_SQUARE_LT_1,                // Wsq (Lorentz factor) is less than 1
    FLUID_NEGATIVE_RHO,                 // primitive density is negative
    FLUID_NEGATIVE_P,                   // primitive variable pressure is negative
    FLUID_V_SQUARE_GT_1,                //  prim, Vsq is GTE 1
    FLUID_MHD_SOLVER_FAILED,            // conserved to primitive MHD solver failed.
    FLUID_MHD_ISENTROPIC_SOLVER_FAILED, // MHD isentropic solver failed.
    FLUID_PRIM_TO_CONS_FAILED,          // overall prim to con var solver failed.
    FLUID_CONS_TO_PRIM_FAILED,          // overall cons to prim solver failed.
    FLUID_CONS_TO_PRIM_RESCALED         // con to prim solver succeeded, but rescaled in the process.
};

/**
 * @brief prints the stack trace for the fluid math code. 
 * 
*/
static inline void print_stacktrace(FILE *out = stderr, unsigned int max_frames = 63)
{
    fprintf(out, "stack trace:\n");

    // storage array for stack trace address data
    void* addrlist[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0) {
	fprintf(out, "  <empty, possibly corrupt>\n");
	return;
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
	char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

	// find parentheses and +address offset surrounding the mangled name:
	// ./module(function+0x15c) [0x8048a6d]
	for (char *p = symbollist[i]; *p; ++p)
	{
	    if (*p == '(')
		begin_name = p;
	    else if (*p == '+')
		begin_offset = p;
	    else if (*p == ')' && begin_offset) {
		end_offset = p;
		break;
	    }
	}

	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';

	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():

	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
					    funcname, &funcnamesize, &status);
	    if (status == 0) {
		funcname = ret; // use possibly realloc()-ed string
		fprintf(out, "  %s : %s+%s\n",
			symbollist[i], funcname, begin_offset);
	    }
	    else {
		// demangling failed. Output function name as a C function with
		// no arguments.
		fprintf(out, "  %s : %s()+%s\n",
			symbollist[i], begin_name, begin_offset);
	    }
	}
	else
	{
	    // couldn't parse the line? print the whole line.
	    fprintf(out, "  %s\n", symbollist[i]);
	}
    }

    free(funcname);
    free(symbollist);
}

/**
 * @brief compute conserved variables form primitive variables. operates on the zip version of the vector. 
 * @param[in] pMesh: underlying mesh data structure. 
 * @param[in] v: primitive variables. (zip version of the vector)
 * @param[out] u : conserved/ evolution vars. (zip version of the vector)
 */
Error prim_to_con(const ot::Mesh *pMesh, double **v, double **u);

/**
 * @brief compute conserved variables form primitive variables for a single grid point
 * @param[out] u : conserved/ evolution vars. 
 * @param[in] v: primitive variables. 
 * @param[in] gd: spacetime metric (down indices)
 * @param[in] sdetg: squaroot of the detg.
 */
Error prim_to_con_pt(double v[], double u[], double gd[3][3], double sdetg);

/***
 * @brief computes the prim vars from conserved variables. (operates on the zipped version. )
 * @param[in] pMesh: underlying mesh data strucutre
 * @param[in/out] u : conserved /evolution variables. 
 * @param[in/out] v : primitive variables. 
 * @note that in the case of the cons to prim solver failure we use interpolation from neighboring points. 
*/
Error con_to_prim(const ot::Mesh *pMesh, double **u, double **v, std::set<unsigned int>& badBounds);

/**
 * @brief compute the primitve vars from constraint variables for a grid point. 
 * @param[in] dupt: conserve variables
 * @param[out] vpt: primitive variable
 * @param[in] gd: metric
 * @param[in] gu: inverse metric
 * @param[in] sdetg: 
 * @param[in] pos: position of the point. 
 * */
Error con_to_prim_pt(double *dupt, double *vpt, double gd[3][3], double gu[3][3], double sdetg, double pos[]);


/**
 * @brief performs 3-vector to 4-vector (velocity) transformation
 * @param[in] v3: 3-vector
 * @param[out] v4: 4-vector
 * @param[in] gd: metric 
 * @param[in] alpha: time like coord 
 * @param[in] beta1: x dir shift vector
 * @param[in] beta2: y dir shift vector
 * @param[in] beta3: z dir shift vector
*/
Error Vvec3_to_Vvec4(double* v3, double* v4, double gd[3][3], double alpha, double Beta[3]);

/**
 * @brief performs 3-vector to 4-vector (velocity) transformation
 * @param[in] v4: 4-vector
 * @param[out] v3: 3-vector
 * @param[in] gd: metric 
 * @param[in] alpha: time like coord 
 * @param[in] beta1: x dir shift vector
 * @param[in] beta2: y dir shift vector
 * @param[in] beta3: z dir shift vector
*/

Error Vvec4_to_Vvec3(double* v4, double* v3, double gd[3][3], double alpha, double Beta[3]);

/**
 * @brief work variables for constraint to primitive solver. 
 */
void c2p_work_vars(double *fpar, double *u, double *v, double gd[3][3], double gu[3][3], double gamma, double *pos);

//@todo: @jacob: Can you please add doxygen documentation tags for the following functions.
int mhd_solver_4(double *u, double *v, double gd[3][3], double gu[3][3], double *fpar);

int mhd_solver_isentropic(double *u, double *v, double *pos, double gd[3][3], double gu[3][3], double *fpar);

int SolvePrimVarsIsentropic(double *qa, double *upt, double *fpar);

int func4_p(double *f, double q, double *fpar);

int func4_p_d(double *f, double *df, double q, double *fpar);

int func4_p_dd(double *f, double *df, double *ddf, double q, double *fpar);

void MathDump(double *fpar);

int func_p_isentropic(double *f, double rho0, double *fpar);

int func_p_d_isentropic(double *f, double *df, double rho0, double *fpar);

} // namespace fluid_math

#endif
