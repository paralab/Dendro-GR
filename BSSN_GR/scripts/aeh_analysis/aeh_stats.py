# compute AEH statistics
import numpy as np
import json
from scipy.special import sph_harm
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 24,
    #"ytick.major.size": 3,
    #"font.family": "Helvetica",
    "lines.linewidth":2.0
})

def sph_harm_real(l, m, theta, phi):
    # in python's sph_harm phi and theta are swapped
    Y = sph_harm(abs(m), l, phi, theta)
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    else:
        Y = Y.real

    return Y 

def load_aeh_data(filename: str):
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
    
    s=",".join(lines)
    s="[" + s + "]"
    return json.loads(s)

def aeh_characteristics(d, ntheta, nphi):
    theta     = np.arccos(np.linspace(-1 , 1, ntheta)) 
    phi       = np.linspace(0 , 2 * np.pi, nphi)
    omega     = np.meshgrid(theta, phi, indexing="ij")
    
    sph_modes = d["lm"]
    h_lm      = np.array(d["coeff"], dtype=np.float64)
    num_sh    = len(sph_modes)
    sph_mat   = np.zeros((ntheta, nphi, num_sh))
    
    for lm_idx, lm in enumerate(sph_modes):
        sph_mat[:, :, lm_idx] = sph_harm_real(lm[0], lm[1], omega[0], omega[1])
    
    r_mat = np.dot(sph_mat, h_lm)
    
    aeh_stat=dict()
    aeh_stat["rmin"]  = np.min(r_mat)
    aeh_stat["rmean"] = np.mean(r_mat)
    aeh_stat["rmax"]  = np.max(r_mat)
    return aeh_stat

def plot_error(aeh_data,fname):
    
    plt.figure(figsize=(30,30),dpi=300)
    
    plt_idx=1
    for d in aeh_data:
        l_max   = d["lm"][-1][0]
        grid    = np.array(d['quad']).reshape((-1, 2))
        
        theta   = np.unique(grid[:,0])
        phi     = np.unique(grid[:,1])
        
        ntheta  = len(theta)
        nphi    = len(phi)
        
        quad    = np.meshgrid(theta, phi, indexing="ij")
        error   = np.abs(np.array(d['expansion']).reshape((ntheta, nphi)))
        
        from matplotlib.colors import LogNorm
        # plt.contour(quad[0], quad[1], error, list(reversed([10**(-i) for i in range(10)])), cmap=plt.cm.jet,norm = LogNorm())
        # plt.xlabel(r"polar angle")
        # plt.ylabel(r"azimuthal angle")
        # plt.colorbar()
        # plt.show()
        
        # from matplotlib.colors import LogNorm
        plt.subplot(2,len(aeh_data)//2 + 1, plt_idx)    
        plt.imshow(np.abs(error), norm=LogNorm(), interpolation='none',extent=[0,np.pi, 0, np.pi * 2], aspect=1)
        plt.xlabel(r"polar angle")
        plt.ylabel(r"azimuthal angle")
        plt.colorbar()
        plt.tight_layout()
        plt.title(r"$l_{max}=%d$"%l_max)
        plt_idx+=1
    plt.suptitle(r"abs(expansion) evaluated on the AH surface")
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    

def plot_convergence(aeh_data,fname):
    plt.figure(figsize=(20,10), dpi=300)
    
    plt.subplot(1,2,1)
    for d in aeh_data:
        lm = d["lm"]
        h  = np.array(d["coeff"])
        plt.semilogy(range(len(h)), np.abs(h), ".-", label=r"l_max=%d"%(lm[-1][0]))
    
    plt.xlabel(r"coefficient id")
    plt.ylabel(r"$|h_{lm}|$")
    plt.grid(visible=True)
    plt.legend()
    
    plt.subplot(1,2,2)
    h_res = np.array(aeh_data[-1]["coeff"])
    rel_error = list()
    lval      = list()
    for d in aeh_data[:-1]:
        lm = d["lm"]
        h  = np.array(d["coeff"])
        
        lval.append(lm[-1][0])
        rel_error.append(np.abs(np.linalg.norm(h)-np.linalg.norm(h_res))/np.linalg.norm(h_res))
        
        
        
    plt.semilogy(np.array(lval, dtype=np.int32), np.array(rel_error) , ".-")
    plt.xlabel(r"$l_{max}$")
    plt.ylabel(r"relative error")
    plt.grid(visible=True)
    plt.savefig(fname)
    plt.close()
    
    
    # plt.show()
    # plt.close()
    
    # plt.subplot(1,3,3)
    
    # lm_res = np.array(aeh_data[-1]["lm"])
    # h_res  = np.array(aeh_data[-1]["coeff"])
    # th     = np.linspace(0,     np.pi , 50)
    # ph     = np.linspace(0, 2 * np.pi , 50)
    # quad   = np.meshgrid(th, ph, indexing="ij")
    
    # aeh_r  = 0
    # for lm_idx, lm in enumerate(lm_res):
    #     aeh_r +=h_res[lm_idx] * sph_harm_real(lm[0], lm[1], quad[0], quad[1])
        
    # plt.contour(quad[0], quad[1], aeh_r)
    
    
    
bh_q4 =[{"rmin":0.08726049192, "rmean": 0.08907037693, "rmax": 0.09093860799}, {"rmin":0.3821575089, "rmean":0.3862445803, "rmax":0.3904143243 }] 

print("\nBH0")
print("AHFinderDirect rmin=%.8E rmean=%.8E rman=%.8E"%(bh_q4[0]["rmin"], bh_q4[0]["rmean"], bh_q4[0]["rmax"]))
runs_folder="../../../build1/q4/dat"    
runs       =[runs_folder+"/dgr_%d_32_32_bh0_aeh.dat"%i for i in range(0,11)]    
aeh_data   =[load_aeh_data(r)[0] for r in runs]
plot_convergence(aeh_data, runs_folder+"/bh0.png")
r_high     = aeh_characteristics(aeh_data[-1],32,32)

for d in aeh_data:
    l_max = d["lm"][-1][0]
    r         = aeh_characteristics(d,32,32)
    rel_error = [abs(r["rmin"]/bh_q4[0]["rmin"]-1) , abs(r["rmean"]/bh_q4[0]["rmean"]-1), abs(r["rmax"]/bh_q4[0]["rmax"]-1)]
    cov_error = [abs(r["rmin"]/r_high["rmin"]-1) , abs(r["rmean"]/r_high["rmean"]-1), abs(r["rmax"]/r_high["rmax"]-1)]
    #print("l max=%d r_min=%.8E r_mean=%.8E r_max=%.8E \t w.r.t AHFinderDirect rel_error(r_min)=%.8E rel_error(r_mean)=%.8E rel_error(r_max)=%.8E"%(l_max, r["rmin"], r["rmean"], r["rmax"], rel_error[0], rel_error[1], rel_error[2]))
    print("%d & %.2E & %.2E & %.2E & %.2E & %.2E & %.2E"%(l_max, rel_error[0], rel_error[1], rel_error[2], cov_error[0], cov_error[1], cov_error[2]))
    

plot_error(aeh_data, runs_folder+"/bh0_expansion.png")


print("\nBH1")
print("AHFinderDirect rmin=%.8E rmean=%.8E rman=%.8E"%(bh_q4[1]["rmin"], bh_q4[1]["rmean"], bh_q4[1]["rmax"]))
runs       =[runs_folder+"/dgr_%d_32_32_bh1_aeh.dat"%i for i in range(0,11)]    
aeh_data   =[load_aeh_data(r)[0] for r in runs]
plot_convergence(aeh_data, runs_folder+"/bh1.png")
r_high     = aeh_characteristics(aeh_data[-1],32,32)
for d in aeh_data:
    l_max = d["lm"][-1][0]
    r         = aeh_characteristics(d,32,32)
    rel_error = [abs(r["rmin"]/bh_q4[1]["rmin"]-1) , abs(r["rmean"]/bh_q4[1]["rmean"]-1), abs(r["rmax"]/bh_q4[1]["rmax"]-1)]
    cov_error = [abs(r["rmin"]/r_high["rmin"]-1) , abs(r["rmean"]/r_high["rmean"]-1), abs(r["rmax"]/r_high["rmax"]-1)]
    #print("l max=%d r_min=%.8E r_mean=%.8E r_max=%.8E \t w.r.t AHFinderDirect rel_error(r_min)=%.8E rel_error(r_mean)=%.8E rel_error(r_max)=%.8E"%(l_max, r["rmin"], r["rmean"], r["rmax"], rel_error[0], rel_error[1], rel_error[2]))
    print("%d & %.2E & %.2E & %.2E & %.2E & %.2E & %.2E"%(l_max, rel_error[0], rel_error[1], rel_error[2], cov_error[0], cov_error[1], cov_error[2]))


plot_error(aeh_data, runs_folder+"/bh1_expansion.png")


