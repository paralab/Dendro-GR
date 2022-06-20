input=inputs[0]
u_chi=input.PointData['U_CHI']
numPoints=input.GetNumberOfPoints()
u_chi_min=min(u_chi)#u_chi[0]
for i in range(numPoints):
    uchi=u_chi[i]
    if abs(uchi-u_chi_min)< 1e-5 :
        minpos=input.GetPoint(i)
        #u_chi_min=uchi
        with open('trajectory.dat', 'a') as ofile:
            ofile.write(str(step)+'\t'+str(minpos[0])+'\t'+str(minpos[1])+'\t'+str(minpos[2])+'\n')
            ofile.close()
        #print "uzmin: ",uchi,min(u_chi)
        #print "minpos: ",minpos
        break;
        #output.append(u_chi_min, 'u_chi_min')
        #output.append(u_chi_min_pos, 'u_chi_min_pos')```
