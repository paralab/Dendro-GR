#include "sc_flops.h"
#include "mpi.h"
#include <stdarg.h>

#ifdef SC_PAPI
#include <sys/types.h>
#include <papi.h>
#endif

void
sc_flops_papi (float *rtime, float *ptime, long long *flpops, float *mflops)
{
#ifdef SC_PAPI
  int                 retval;

  retval = PAPI_flops (rtime, ptime, flpops, mflops);
  SC_CHECK_ABORT (retval == PAPI_OK, "Papi not happy");
#else
  *rtime = *ptime = *mflops = 0.;
  *flpops = 0;
#endif
}

void
sc_flops_start (sc_flopinfo_t * fi)
{
  float               rtime, ptime, mflops;
  long long           flpops;

  fi->seconds = MPI_Wtime ();
  sc_flops_papi (&rtime, &ptime, &flpops, &mflops);     /* ignore results */

  fi->cwtime = 0.;
  fi->crtime = fi->cptime = 0.;
  fi->cflpops = 0;

  fi->iwtime = 0.;
  fi->irtime = fi->iptime = fi->mflops = 0.;
  fi->iflpops = 0;
}

void
sc_flops_count (sc_flopinfo_t * fi)
{
  double              seconds;
  float               rtime, ptime;
  long long           flpops;

  seconds = MPI_Wtime ();
  sc_flops_papi (&rtime, &ptime, &flpops, &fi->mflops);

  fi->iwtime = seconds - fi->seconds;
  fi->cwtime += fi->iwtime;

  fi->iptime = ptime - fi->cptime;
  fi->cptime = ptime;

  fi->iflpops = flpops - fi->cflpops;
  fi->cflpops = flpops;

#ifdef SC_PAPI
  fi->irtime = rtime - fi->crtime;
  fi->crtime = rtime;
#else
  fi->irtime = (float) fi->iwtime;
  fi->crtime = (float) fi->cwtime;
#endif
  fi->seconds = seconds;
}

void
sc_flops_snap (sc_flopinfo_t * fi, sc_flopinfo_t * snapshot)
{
  sc_flops_count (fi);
  *snapshot = *fi;
}

void
sc_flops_shot (sc_flopinfo_t * fi, sc_flopinfo_t * snapshot)
{
  sc_flops_shotv (fi, snapshot, NULL);
}

void
sc_flops_shotv (sc_flopinfo_t * fi, ...)
{
  sc_flopinfo_t      *snapshot;
  va_list             ap;

  sc_flops_count (fi);

  va_start (ap, fi);
  for (; (snapshot = va_arg (ap, sc_flopinfo_t *)) != NULL;) {
    snapshot->iwtime = fi->cwtime - snapshot->cwtime;
    snapshot->irtime = fi->crtime - snapshot->crtime;
    snapshot->iptime = fi->cptime - snapshot->cptime;
    snapshot->iflpops = fi->cflpops - snapshot->cflpops;
    snapshot->mflops =
      (float) ((double) snapshot->iflpops / 1.e6 / snapshot->irtime);

    snapshot->seconds = fi->seconds;
    snapshot->cwtime = fi->cwtime;
    snapshot->crtime = fi->crtime;
    snapshot->cptime = fi->cptime;
    snapshot->cflpops = fi->cflpops;
  }
  va_end (ap);
}
