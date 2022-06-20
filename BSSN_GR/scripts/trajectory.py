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
    hadName=sys.argv[2]

    xmin=-8
    xmax=8

    ymin=-8
    ymax=8

    tstep=list()
    tstime=list()
    bh1_x=list()
    bh1_y=list()
    bh1_z=list()

    bh2_x=list()
    bh2_y=list()
    bh2_z=list()


    hX=list()
    hY=list()
    hT=list()

    with open(hadName) as tsv:
        reader = csv.reader(tsv)
        for line in csv.reader(tsv, delimiter='\t'):
            #line = list(filter(None, line))
            #line=str(line).split()
            line = [item for item in filter(None, line[0].split(' '))]
            #print(line)
            hT.append(float(line[0]))
            hX.append(float(line[2]))
            hY.append(float(line[3]))


    with open(fName) as tsv:
        reader = csv.reader(tsv)
        next(reader,None)
        for line in csv.reader(tsv, delimiter='\t'):
            line = list(filter(None, line))
            #print(line)
            tstep.append(float(line[0]))
            tstime.append(float(line[0])*(0.25*(400.0/(2**14))))
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


    hX=np.array(hX)
    hY=np.array(hY)
    hT=np.array(hT)

    
    fig = plt.figure();
    plt.subplot(2,2,1)
    #plt.scatter(bh1_x, bh1_y,label='BH1',s=1,'r-')
    #plt.scatter(bh2_x, bh2_y,label='BH2',s=1,'b-')
    plt.plot(bh1_x, bh1_y,label='BH1',linestyle='-')
    plt.plot(bh2_x, bh2_y,label='BH2',linestyle='-')
    plt.legend(loc='upper left')
    plt.axis('equal')
    plt.title('x,y plot')
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])
    #plt.show()

    plt.subplot(2,2,2)
    #plt.scatter(tstep,bh1_z,label='BH_1_z',s=1)
    #plt.scatter(tstep,bh2_z,label='BH_2_z',s=1)
    plt.plot(tstime,bh1_z,label='BH_1_z',linestyle='-')
    plt.plot(tstime,bh2_z,label='BH_2_z',linestyle='-')
    plt.title('z axis')
 
 
    plt.subplot(2,2,3)
    distHad=np.sqrt((2*hX)**2 + (2*hY)**2)
    distDendro=np.sqrt((bh1_x-bh2_x)**2 + (bh1_y-bh2_y)**2)
    plt.plot(hT,distHad,label='had',linestyle='-')
    plt.plot(tstime,distDendro,label='dendro',linestyle='-')
    plt.legend(loc='upper right')
    plt.title('distance')
    axes = plt.gca()
    #axes.set_xlim([xmin,xmax])
    axes.set_ylim([0,9])

    plt.subplot(2,2,4)
    plt.plot(bh1_x, bh1_y,label='dendro',linestyle='-')
    plt.plot(hX, hY,label='HAD',linestyle='-')
    plt.legend(loc='upper right')
    plt.axis('equal')
    plt.title('trajectory')
    axes = plt.gca()
    axes.set_xlim([xmin,xmax])
    axes.set_ylim([ymin,ymax])

    
 
    fig.savefig('trajectory.png')
    plt.close(fig)



    
if __name__ == "__main__":
    main()

