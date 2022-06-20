# author: Milinda Fernando
# script to plot the trajectory of the black holes with different dendro parameters. 
# School of Computing, University of Utah
# 09/13/2019

import argparse
import resource
import gc
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
import csv as csv

def main():
    # print command line arguments
    if(len(sys.argv)==0):
        print("Error: bh location file not specified. ")
        sys.exit(0)
    
    fName1=sys.argv[1]
    fName2=sys.argv[2]

    xmin=-8
    xmax=8

    ymin=-8
    ymax=8

    tstep1=list()
    tstime1=list()
    
    tstep2=list()
    tstime2=list()

    f1_bh1_x=list()
    f1_bh1_y=list()
    f1_bh1_z=list()

    f1_bh2_x=list()
    f1_bh2_y=list()
    f1_bh2_z=list()

    f2_bh1_x=list()
    f2_bh1_y=list()
    f2_bh1_z=list()

    f2_bh2_x=list()
    f2_bh2_y=list()
    f2_bh2_z=list()


    with open(fName1) as tsv:
        reader = csv.reader(tsv)
        next(reader,None)
        for line in csv.reader(tsv, delimiter='\t'):
            line = list(filter(None, line))
            #print(line)
            tstep1.append(float(line[0]))
            tstime1.append(float(line[0])*(0.25*(400.0/(2**14))))
            f1_bh1_x.append(float(line[1]))
            f1_bh1_y.append(float(line[2]))
            f1_bh1_z.append(float(line[3]))

            f1_bh2_x.append(float(line[4]))
            f1_bh2_y.append(float(line[5]))
            f1_bh2_z.append(float(line[6]))

    f1_bh1_x=np.array(f1_bh1_x)
    f1_bh1_y=np.array(f1_bh1_y)
    f1_bh1_z=np.array(f1_bh1_z)
    
    f1_bh2_x=np.array(f1_bh2_x)
    f1_bh2_y=np.array(f1_bh2_y)
    f1_bh2_z=np.array(f1_bh2_z)

    with open(fName2) as tsv:
        reader = csv.reader(tsv)
        next(reader,None)
        for line in csv.reader(tsv, delimiter='\t'):
            line = list(filter(None, line))
            #print(line)
            tstep2.append(float(line[0]))
            tstime2.append(float(line[0])*(0.25*(400.0/(2**14))))
            f2_bh1_x.append(float(line[1]))
            f2_bh1_y.append(float(line[2]))
            f2_bh1_z.append(float(line[3]))

            f2_bh2_x.append(float(line[4]))
            f2_bh2_y.append(float(line[5]))
            f2_bh2_z.append(float(line[6]))

    f2_bh1_x=np.array(f2_bh1_x)
    f2_bh1_y=np.array(f2_bh1_y)
    f2_bh1_z=np.array(f2_bh1_z)
    
    f2_bh2_x=np.array(f2_bh2_x)
    f2_bh2_y=np.array(f2_bh2_y)
    f2_bh2_z=np.array(f2_bh2_z)

    fig = plt.figure()

    plt.subplot(2,1,1)
    distDendro1=np.sqrt((f1_bh1_x-f1_bh2_x)**2 + (f1_bh1_y-f1_bh2_y)**2)
    distDendro2=np.sqrt((f2_bh1_x-f2_bh2_x)**2 + (f2_bh1_y-f2_bh2_y)**2)

    plt.plot(tstime1,distDendro1,label='dendro1',linestyle='-')
    plt.plot(tstime2,distDendro2,label='dendro2',linestyle='-')
    plt.legend(loc='upper right')
    plt.title('distance')
    axes = plt.gca()
    #axes.set_xlim([xmin,xmax])
    axes.set_ylim([0,9])

    fig.savefig('ctrajectory.png')
    plt.close(fig)

if __name__ == "__main__":
    main()

















