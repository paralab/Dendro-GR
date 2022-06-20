# author: Milinda Fernando
# script to plot the trajectory of the black holes
# School of Computing, University of Utah
# 01/17/2019
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
    
    fName=sys.argv[1]
    xmin=-8
    xmax=8

    ymin=-8
    ymax=8

    tstep=list()
    bh1_x=list()
    bh1_y=list()
    bh1_z=list()

    bh2_x=list()
    bh2_y=list()
    bh2_z=list()

    with open(fName) as tsv:
        reader = csv.reader(tsv)
        next(reader,None)
        for line in csv.reader(tsv, delimiter='\t'):
            line = list(filter(None, line))
            #print(line)
            tstep.append(float(line[0]))
            bh1_x.append(float(line[1]))
            bh1_y.append(float(line[2]))
            bh1_z.append(float(line[3]))

            bh2_x.append(float(line[4]))
            bh2_y.append(float(line[5]))
            bh2_z.append(float(line[6]))

    bh1_x=np.array(bh1_x)
    bh1_y=np.array(bh1_y)
    bh1_z=np.array(bh1_z)

    bh2_x=np.array(bh2_x)
    bh2_y=np.array(bh2_y)
    bh2_z=np.array(bh2_z)
    
    plt.scatter(bh1_x, bh1_y,label='BH1')
    plt.scatter(bh2_x, bh2_y,label='BH2')
    plt.legend(loc='upper left')
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    plt.show()

    
if __name__ == "__main__":
    main()

