import math as m
import cmath as c
import numpy as np

import gw

gw_dendro = gw.GWExtract(5,131,[2,3,4],2,[])

def dendro_swsh(Y, theta, phi):
    Y[0] = gw_dendro.swsh(-2, 2, -2, theta, phi)
    Y[1] = gw_dendro.swsh(-2, 2, -1, theta, phi)
    Y[2] = gw_dendro.swsh(-2, 2,  0, theta, phi)
    Y[3] = gw_dendro.swsh(-2, 2,  1, theta, phi)
    Y[4] = gw_dendro.swsh(-2, 2,  2, theta, phi)


def nrswsh(Y, theta, phi):
    Y[0] = np.sqrt(5.0/(64.0 *np.pi))*((1.0 - np.cos(theta))**2)*np.exp(-2.0*(0+1.0j*phi))
    Y[1] = np.sqrt(5.0/(16.0 *np.pi))*np.sin(theta)*(1.0 - np.cos(theta))*np.exp(-(0+1.0j*phi))
    Y[2] = np.sqrt(15.0/(32.0*np.pi))*np.sin(theta)**2
    Y[3] = np.sqrt(5.0/(16.0 *np.pi))*np.sin(theta)*(1.0 + np.cos(theta))*np.exp((0+1.0j*phi))
    Y[4] = np.sqrt(5.0/(64.0 *np.pi))*((1.0 + np.cos(theta))**2)*np.exp(2.0*(0+1.0j*phi))



theta = np.array(gw_dendro.lebedev_theata)
phi   = np.array(gw_dendro.lebedev_phi)

Y_dendro_list=[0,0,0,0,0]
Y_nr_list =[0,0,0,0,0]

dendro_swsh(Y_dendro_list,theta,phi)
dendro_swsh(Y_nr_list,theta,phi)

Y_diff = np.abs(np.abs(Y_dendro_list[1])-np.abs(Y_nr_list[1]))
print(np.max(Y_diff))
