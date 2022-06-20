'''
@author: Milinda Fernando
@brief: Computes the strain of GWs. 
        Based on getStarin.py shared by David. 

For more details on GW extraction see[1].

[1]. Bishop, N. T., & Reisswig, C. (2014). The gravitational wave strain in the characteristic
formalism of numerical relativity. General Relativity and Gravitation, 46(1), 1643.
'''

import numpy as np
import matplotlib.pyplot as pl
from pylab import *
import scipy as scipy
from scipy import constants
from scipy.fftpack import fft
from scipy.interpolate import interp1d
from scipy.integrate   import cumtrapz
import pandas as pd
import sys as sys
from scipy.interpolate import make_interp_spline, BSpline

'''
compute the propergation time. based in schwarchild metric
'''
def compute_tstar(t,r,M=1.0,c=1.0):
   '''
   compute t^*(t,s) =\int_{0}^{t^*} \frac{\sqrt{-g^{ss}/ g^{tt}}}{1-2M/r} dt^\prime =r -2M ln (r/2M -1)
   '''
   ts = t  - ( r + 2*M*log((r/2*M)-1) )
   return ts   


'''
compute t from the t* and specific r value. 
'''
def compute_t(ts,r,M=1.0,c=1.0):
   if(r/2*M<=1):
      print("r/2M<1 hence t is not physical")
      t=0
   else:
      t = ts + ( r + 2*M*log((r/2*M)-1) )

   return t



'''
Read the psi4 lm values from file. 
@param [in] fprefix: file prefix. 
@param [in] l : l mode
@param [in] m : m mode
@param [in] k : number of different radaii
@return [t,psi4]
'''
def read_data(fprefix,l,m,k):
   
   fname= fprefix + "_l" + str(l) + "_m" + str(abs(m))+ ".dat"
   df=pd.read_csv(fname,sep='\t')
   df.columns = df.columns.str.replace(' ', '')
   
   t=np.array(df["t"])
   psi4=list()
   for r in range(0,k):
      psi4_r=list()
      for item in df["r"+str(r)]:
         [x,y]=str(item).strip('()').split(',')
         psi4_r.append(float(x) + 1j* float(y))
      psi4.append(np.array(psi4_r))

   return [t,psi4]


'''
compute psi4 at infinity for constant t* values
'''
def get_psi4_infinity(t,psi4,r):
   r_min=min(r)
   r_max=max(r)

   numr=len(r)

   t_min=min(t)
   t_max=max(t)

   psi4T_real = np.zeros((len(r),len(t)))
   psi4T_img  = np.zeros((len(r),len(t)))
   
   for i in range(0,len(r)):
      psi4T_real[i,:] = np.real(psi4[i])
      psi4T_img [i,:] = np.imag(psi4[i])

   # interpolated psi4 with t and r
   fpsi4_real = scipy.interpolate.interp2d(t,r,psi4T_real,kind='cubic')
   fpsi4_img = scipy.interpolate.interp2d(t,r,psi4T_img,kind='cubic')

   s = 1.0/np.linspace(1.0/r_max, 1.0/r_min , numr )
   s=sort(s)
   
   #print(s)
   
   ts_max = compute_tstar(t_max,r_min)
   ts = np.arange(0,ts_max,0.5)

   tt = np.zeros(numr)

   psi4_ts = np.zeros((len(ts),numr) , dtype=np.complex)
   psi4_ts_inf = np.zeros(len(ts) , dtype=np.complex)

   for i in range(0,len(ts)):
      for j in range(0,numr):
         tt[j] = compute_t(ts[i],s[j])
         psi4_ts[i,j] = fpsi4_real(tt[j],s[j]) + 1.0j * fpsi4_img(tt[j],s[j])

   #print(np.size(psi4_ts,0))
   #print(np.size(psi4_ts,1))
   #print(psi4_ts)   

   V = np.zeros((numr,numr))
   rhs_abs = np.zeros((numr,1))
   rhs_arg = np.zeros((numr,1))

   coeff_abs =np.zeros((numr,1))
   coeff_arg =np.zeros((numr,1))

   for i in range(0,numr):
      for j in range(0,numr):
         V[i,j] = 1.0/(s[i]**j)

   invV=np.linalg.inv(V)

   for i in range(0,len(ts)):
      for j in range(0,numr):
         rhs_abs[j,0]=np.abs(psi4_ts[i,j])
         rhs_arg[j,0]=np.angle(psi4_ts[i,j])

      coeff_abs = np.matmul(invV,rhs_abs)
      coeff_arg = np.matmul(invV,rhs_arg)

      psi4_ts_inf[i] = coeff_abs[0,0] * (cos(coeff_arg[0,0]) + 1j * sin (coeff_arg[0,0]))

   #print(psi4_ts_inf)
   return [ts,psi4_ts_inf]

##################################################
##     strain --via FFT->-1/omega^2->iFFT       ##
##################################################
def strain(ts, psi4_inf, cutoff):
   
   # Compute the FFT:
   fftU   = fft(psi4_inf)
   tstep  = ts[1]-ts[0]
   N      = len(ts)
   
   # Compute the frequencies of the FFT:
   fftU_f = fftfreq(N,d=tstep)
   # Kill the DC offset:
   fftU[0]= 0.0
   # Change the low-end frequencies
   freqs=fftU_f
   # Make sure zero frequency does not cause error
   freqs[0]=1.0
   for i in range(len(fftU_f)):
      if (abs(fftU_f[i]) < cutoff):
         freqs[i] = cutoff*sign(fftU_f[i])
      else:
         freqs[i] = fftU_f[i]
   
   # Do the integration and inverse FFT:
   h = -ifft(fftU/((2.0*pi*freqs)**2))
   return h

def plt_psi4(fname,lmodes,fprefix):
   df=pd.read_csv(fname,sep='\t')
   df.columns = df.columns.str.replace(' ', '')
   t=np.array(df["t"])
   
   for l in lmodes:

      fig = plt.figure()
      for m in range(-l,l+1):
         if m >=0: 
            colabs="abs(l_%dm_p_%d)"%(l,abs(m))
            colarg="arg(l_%dm_p_%d)"%(l,abs(m))
         else:
            colabs="abs(l_%dm_m_%d)"%(l,abs(m))
            colarg="arg(l_%dm_m_%d)"%(l,abs(m))

         r=np.array(df[colabs])
         theta=np.array(df[colarg])

         psi4_real=r*cos(theta)
         psi4_img=r*sin(theta)

         fR = interp1d(t, psi4_real, kind='cubic')
         fI = interp1d(t, psi4_img, kind='cubic')
         tnew=np.linspace(t.min(),t.max(),num=10000,endpoint=True)

         subplot(2*l+1,2,2*(l+m)+1)
         plt.plot(tnew,fR(tnew))
         plt.title("real l=%d,m=%d"%(l,m))
         plt.gcf().set_size_inches(50, 18)
         
         subplot(2*l+1,2,2*(l+m)+2)
         plt.plot(tnew,fI(tnew))
         plt.title("imag l=%d,m=%d"%(l,m))
         plt.gcf().set_size_inches(50, 18)

      fig.savefig("%s_l%d.png"%(fprefix,l))
      plt.close(fig)

   

def plt_strain(ts,hstrain,l,m,fprefix="GW_strain"):
   fig = plt.figure()
   
   subplot(2,1,1)
   plt.plot(ts,np.real(hstrain))
   plt.title("h_+ l=%d,m=%d"%(l,m))
   plt.gcf().set_size_inches(50, 9)

   subplot(2,1,2)
   plt.plot(ts,np.imag(hstrain))
   plt.title("h_x l=%d,m=%d"%(l,m))
   plt.gcf().set_size_inches(50, 9)

   fig.savefig("%s_.png"%(fprefix))
   plt.close(fig)



def main():
   # print command line arguments
   if(len(sys.argv) == 0):
      print("Error: psi4 data file is not provided ")
      sys.exit(0)

   fprefix=sys.argv[1]
   #plt_psi4(fname,[2,3,4],"r1")
   #plt_strain(fname,[2,3,4],"r1_strain")
   
   l=2
   m=-2
   
   [t,psi4] = read_data(fprefix,l,m+l,5)
   [ts,psi4_inf]=get_psi4_infinity(t,psi4,[50,70,90,110,130])
   
   hstrain= strain(ts,psi4_inf,2)
   
   plt_strain(ts,hstrain,l,m,"GW_l%d_m%d" %(l,m))
   
if __name__ == "__main__":
    main()





