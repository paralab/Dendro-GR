#ifndef SC_FLOPS_H
#define SC_FLOPS_H

#include "stddef.h"

typedef struct sc_flopinfo
{
  double              seconds;  /* current time from MPI_Wtime */

  /* these variables measure onward from from sc_flops_start */
  double              cwtime;   /* cumulative wall time */
  float               crtime;   /* cumulative real time */
  float               cptime;   /* cumulative process time */
  long long           cflpops;  /* cumulative floating point operations */

  /* measure since sc_flops_start or the previous sc_flops_count */
  double              iwtime;   /* interval wall time */
  float               irtime;   /* interval real time */
  float               iptime;   /* interval process time */
  long long           iflpops;  /* interval floating point operations */
  float               mflops;   /* MFlop/s rate in this interval */

  /* without SC_PAPI only seconds, ?wtime and ?rtime are meaningful */
}
sc_flopinfo_t;

/**
 * Calls PAPI_flops.  Aborts on PAPI error.
 * The first call sets up the performance counters.
 * Subsequent calls return cumulative real and process times,
 * cumulative floating point operations and the flop rate since the last call.
 */
void                sc_flops_papi (float *rtime, float *ptime,
                                   long long *flpops, float *mflops);

/**
 * Prepare sc_flopinfo_t structure and start flop counters.
 * Must only be called once during the program run.
 * This function calls sc_flops_papi.
 *
 * \param [out] fi  Members will be initialized.
 */
void                sc_flops_start (sc_flopinfo_t * fi);

/**
 * Update sc_flopinfo_t structure with current measurement.
 * Must only be called after sc_flops_start.
 * Can be called any number of times.
 * This function calls sc_flops_papi.
 *
 * \param [in,out] fi   Members will be updated.
 */
void                sc_flops_count (sc_flopinfo_t * fi);

/**
 * Call sc_flops_count (fi) and copies fi into snapshot.
 *
 * \param [in,out] fi       Members will be updated.
 * \param [out] snapshot    On output is a copy of fi.
 */
void                sc_flops_snap (sc_flopinfo_t * fi,
                                   sc_flopinfo_t * snapshot);

/**
 * Call sc_flops_count (fi) and override snapshot interval timings
 * with the differences since the previous call to sc_flops_snap.
 * The interval mflop rate is computed by iflpops / 1e6 / irtime.
 * The cumulative timings in snapshot are copied form fi.
 *
 * \param [in,out] fi       Members will be updated.
 * \param [in,out] snapshot Interval timings measured since sc_flops_snap.
 */
void                sc_flops_shot (sc_flopinfo_t * fi,
                                   sc_flopinfo_t * snapshot);

/**
 * Call sc_flops_count (fi) and work on all arguments in the list
 * of type sc_flopinfo_t * as in sc_flops_shot.  Last argument must be NULL.
 */
void                sc_flops_shotv (sc_flopinfo_t * fi, ...);


#endif /* !SC_FLOPS_H */
