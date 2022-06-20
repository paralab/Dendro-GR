import sys
import numpy as np
import lalsimutils
import lal
import lalsimulation as lalsim

import geometric_to_MKS
import re
import argparse

SEG_LENGTH=4.0  # length of segment in seconds
TIME_DIFF_TOL = 1.0e-9

def clean_time_series(time, data, assume_uniform_grid=False):
    """
    PURPOSE: After a simulation restart, time series data already
    output may be outputed again. This routine looks for
    non-monotinicity in the time series.
    Inputs:
        time: a np.ndarray of times
        data: a np.ndarray of data values corresponding to the time
        assume_uniform_grid: Boolean that specifies whether or not the time
        series should be assumed to be uniform. When set to true, the correct
        (cleaned) time series is assumed to be time[i]= time[0] + i *
        (time[1]-time[0]). Any deviation from this is assumed to be due to
        truncation error in the time array and "fixed".
     """
    
    if type(time) != np.ndarray:                                    
        sys.stderr.write("clean_time_series: requires an ndarray\n")
        sys.exit(-1)

    if type(data) != np.ndarray:                                    
        sys.stderr.write("clean_time_series: requires an ndarray\n")
        sys.exit(-1)

    ntime = []
    ndata = []
    ntime.append(time[0])
    ndata.append(data[0])
    ntime.append(time[1])
    ndata.append(data[1])

    expected_dt = time[1] - time[0]
    last_time = time[1]
    tcounter = 1
    if (assume_uniform_grid):
        for i in range(2, len(time)):
            if time[i] - last_time > expected_dt - expected_dt/5:
                tcounter += 1
                last_time = expected_dt * tcounter + time[0]
                ntime.append(expected_dt * tcounter + time[0])
                ndata.append(data[i])
                if np.abs(time[i] - time[0] - tcounter * expected_dt) > expected_dt/5.:
                    sys.stderr.write("Non uniform grid. Fatal Error\n")
                    sys.exit(-1)

            else:
                sys.stderr.write("Skipped entry : {}\n".format(time[i]))
    else:
        for i in range(2, len(time)):
            if time[i] - last_time > TIME_DIFF_TOL:
                last_time = time[i]
                ntime.append(time[i])
                ndata.append(data[i])
            else:
                sys.stderr.write("Skipped entry : {}\n".format(time[i]))
    ntime = np.array(ntime)
    ndata = np.array(ndata)

    if len(ndata) != len(data):
        sys.stderr.write("Time series cleaned: redundant data discarded\n")
        sys.stderr.write("This is not a fatal error\n")

    return ntime, ndata

def nonuniform_to_uniform_grid(time, data):
    """
    PURPOSE: This routine takes a monotonic time series and interpolated it
    onto a unform grid of the form new_time[i] = time[0] +  i * dt,
    where dt is the minimun time difference between neighboring data points
    in the "time" array.

    new_time[-1] lies inside (time[-1]-dt, time[-1]]

    Inputs:
        time: a np.ndarray of times
        data: a np.ndarray of data values corresponding to the time
    """

    if type(time) != np.ndarray:                                    
        sys.stderr.write("clean_time_series: requires an ndarray\n")
        sys.exit(-1)

    if type(data) != np.ndarray:                                    
        sys.stderr.write("clean_time_series: requires an ndarray\n")
        sys.exit(-1)

    mindt = np.min(time[1:]-time[0:-1])
    npts = int((time[-1] - time[0])/mindt)
    ntime = np.linspace(time[0], time[0] + mindt*npts, npts+1)
    ndata = np.interp(ntime, time, data)
    return ntime, ndata

def nextPow2(length):        
    """                                                  
    Find next power of 2 <= length of data                                                                              
    """                                                    
    return int(2**np.ceil(np.log2(length)))

def str2z(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return complex(float(sr), float(si))

def str2r(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return float(sr)

def str2i(s):
    s = s.decode("utf-8")
    ss = re.sub('[\(\)]','',s)
    sr, si = ss.split(',')
    return float(si)

def overlap(wave1, wave2, mass):
    """
    Purpose: Calculate the overlap of two waveforms. Both waveforms are assumed
    to be psi_4 rather than the strain.
    INPUTS:
        wave1, wave2: Dictionaries with entries "time", and "psi4" corresponding to
            a uniform time grid and the complex psi_4, respectively.
        mass: the total mass in units of solar masses
    """

    # The mass normalization only affects the time array (any multiplicative
    # change to psi_4 is canceled and can therefore be ignored

    wave1_nt = geometric_to_MKS.GeometricTime_To_MKS_Time(wave1["time"],
            mass*geometric_to_MKS.msun)
    wave2_nt = geometric_to_MKS.GeometricTime_To_MKS_Time(wave2["time"],
            mass*geometric_to_MKS.msun)
  

    # both waveforms are interpolated onto a standard grid
    t_final = np.min((wave1_nt[-1], wave2_nt[-1]))
    dt = np.min((wave1_nt[1]-wave1_nt[0], wave2_nt[1]-wave2_nt[0]))
    npts = int(t_final / dt)
  
    tgrid = np.linspace(0.0, npts * dt, npts+1)
    wave1_psi4_grid = np.interp(tgrid, wave1_nt, wave1["psi4"])
    wave2_psi4_grid = np.interp(tgrid, wave2_nt, wave2["psi4"])
  
    Psi4T1 = lal.CreateCOMPLEX16TimeSeries("Complex overlap",lal.LIGOTimeGPS(0.), 0., dt, lal.DimensionlessUnit, npts+1)
    Psi4T2 = lal.CreateCOMPLEX16TimeSeries("Complex overlap",lal.LIGOTimeGPS(0.), 0., dt, lal.DimensionlessUnit, npts+1)
  
    Psi4T1.data.data = wave1_psi4_grid
    Psi4T2.data.data = wave2_psi4_grid
  
    # The final time is set by the desired segment length and the restriction that
    # the number of points must be a power of 2. We enlarge the grids by adding
    # zeros
    nsize = int(SEG_LENGTH / dt)
    pow2npts = nextPow2(nsize)
    Psi4T1 = lal.ResizeCOMPLEX16TimeSeries(Psi4T1, 0, pow2npts)
    Psi4T2 = lal.ResizeCOMPLEX16TimeSeries(Psi4T2, 0, pow2npts)


    Psi4F1 = lalsimutils.DataFourier(Psi4T1)
    Psi4F2 = lalsimutils.DataFourier(Psi4T2)
  
  
    Tlen = pow2npts * dt
    psd=lalsim.SimNoisePSDaLIGOZeroDetHighPower
    IP = lalsimutils.CreateCompatibleComplexOverlap(Psi4F1,
          psd=psd,fLow=20,fMax=2000,analyticPSD_Q=True,interpolate_max=False,
          waveform_is_psi4=True)
  
    rho_1 = IP.norm(Psi4F1)
    rho_2 = IP.norm(Psi4F2)
    inner_12 = IP.ip(Psi4F1,Psi4F2)/rho_1/rho_2
    return inner_12

def read_dendro_waveform(wave_dir, obs, l, m):
    # Load DendroGR data
    # Can't use str2z directly due to incompatibility with Transpose
    dendro_re  = np.genfromtxt(wave_dir+"/bssn_prof_GW_l{0}_m{1}.dat".format(L, M),
            converters={2: str2r, 3: str2r, 4: str2r, 5: str2r, 6 : str2r, 7 :
            str2r}, skip_header=1).T
    dendro_im  = np.genfromtxt(wave_dir+"/bssn_prof_GW_l{0}_m{1}.dat".format(L, M),
            converters={2: str2i, 3: str2i, 4: str2i, 5: str2i, 6 : str2i, 7 :
            str2i}, skip_header=1).T

    de_time = dendro_re[1]
    de_psi4 = dendro_re[obs + 2]  + 1j * dendro_im[obs +2]

    time, psi4 = nonuniform_to_uniform_grid(*clean_time_series(de_time, de_psi4))
    return {"time": time, "psi4": psi4}

def read_lazev_waveform(wave_dir, obs, l, m):
    # Load LazEv data        
    le_time, le_real, le_imag, _, _ =\
            np.loadtxt(wave_dir+"/pk_rad_obs_{0}_psi4_{1}_{2}.tl".format(OBS, L, M)).T
    le_psi4 = le_real + 1j * le_imag

    time, psi4 = nonuniform_to_uniform_grid(*clean_time_series(le_time, le_psi4))
    return {"time": time, "psi4": psi4}



parser = argparse.ArgumentParser()
parser.add_argument("--dir1", dest='dir1',\
        help="Directory containing first waveform output files", type=str)
parser.add_argument("--waveform1-code", dest='wave1_code',\
        help="Code that generated waveform 1", type=str)
parser.add_argument("--waveform2-code", dest='wave2_code',\
        help="Code that generated waveform 2", type=str)

parser.add_argument("--dir2", dest='dir2',\
        help="Directory containing second waveform-output files", type=str)
parser.add_argument("--l-mode", dest='lmode', help="l-mode", type=int)
parser.add_argument("--m-mode", dest='mmode', help="m-mode", type=int)
parser.add_argument("--observer", help="observer number", type=int)

parser.set_defaults(dir1=".")
parser.set_defaults(dir2=".")
parser.set_defaults(lmode=2)
parser.set_defaults(mmode=2)
parser.set_defaults(observer=0)
parser.set_defaults(wave1_code="lazev")
parser.set_defaults(wave2_code="dendro")

args = parser.parse_args()

OBS=args.observer
L=args.lmode
M=args.mmode

if args.wave1_code == "lazev":
    wave1 = read_lazev_waveform(args.dir1, OBS, L, M)
elif args.wave1_code == "dendro":
    wave1 = read_dendro_waveform(args.dir1, OBS, L, M)
else:
    sys.stderr.write("Error. Unknown waveform1-code: {}\n".format(args.wave1_code))
    sys.exit(-1)

if args.wave2_code == "lazev":
    wave2 = read_lazev_waveform(args.dir2, OBS, L, M)
elif args.wave2_code == "dendro":
    wave2 = read_dendro_waveform(args.dir2, OBS, L, M)
else:
    sys.stderr.write("Error. Unknown waveform2-code: {}\n".format(args.wave2_code))
    sys.exit(-1)



print ("# Comparison of data in {} and {} for l={}, m={}, obs={}".format(args.dir1, args.dir2, L, M, OBS))

for mass in range(10, 205, 5):
    print (mass, overlap(wave1, wave2, mass))
