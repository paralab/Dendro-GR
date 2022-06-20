#!/usr/bin/env python
#####
### CS6635 - Scientific Visualization Project. VTK based visualization framework for computational relativity. 
### Author: Milinda Fernando, Max Carlson
### School of Computing, University of Utah.
### Date: 04/02/2018
####

'''

All the filters, and filter helper functions. 

'''

import vtk as vtk
from mpi4py import MPI
import numpy as np

'''
compute the slice of the octree. 
'''

def SliceFilter(source,genScalars=True,genTriangles=False,origin=[2048,2048,2048],normal=[0,0,1]):

    plane=vtk.vtkPlane()
    plane.SetOrigin(origin[0],origin[1],origin[2])
    plane.SetNormal(normal[0],normal[1],normal[2])

    #create cutter
    cutter=vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(source.GetOutputPort())
    # genScalars on computes the scalar values at the plane. 
    
    if genScalars:  
        cutter.GenerateCutScalarsOn()
    else:
        cutter.GenerateCutScalarsOff()
    
    #genTriangles on compute the triangular mesh on the slice
    
    '''if genTriangles:
        cutter.GenerateTrianglesOn()
    else:
        cutter.GenerateTrianglesOff()
    '''
    cutter.Update()
    
    return cutter


'''
color by a specific variable or a vector. 
'''
def ComputeScalarRange(source,varName):
    pointData=source.GetOutput().GetPointData()
    data=pointData.GetArray(varName)
    data_range=data.GetRange()
    return data_range


'''
warp by scalar
'''
def WarpByScalar(source,varName,scaleFactor=10.0):
    warpByScalar=vtk.vtkWarpScalar()
    warpByScalar.SetInputConnection(source.GetOutputPort())
    warpByScalar.SetScaleFactor(scaleFactor); #use the scalars themselves
    warpByScalar.UseNormalOn();
    warpByScalar.SetNormal(0, 0, 1);
    warpByScalar.Update()
    return warpByScalar


'''
Generate Diverging color scheme with RGB1 and RGB2.

'''

def GenerateDivergingColorMap(RGB1=[0,0,1],RGB2=[1,0,0],numColors=256):
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(numColors)
    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToDiverging()
    ctf.AddRGBPoint(0,float(RGB1[0]),float(RGB1[1]),float(RGB1[2]))
    ctf.AddRGBPoint(1,float(RGB2[0]),float(RGB2[1]),float(RGB2[2]))
    #ctf.AddRGBPoint(0.0, 0, 0, 1.0)
    #ctf.AddRGBPoint(1.0, 1.0, 0, 0)
    for ii,ss in enumerate([float(xx)/float(numColors) for xx in range(numColors)]):
        cc = ctf.GetColor(ss)
        lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
    return lut



'''
Extract the BH locations based on the shortest L2 distance to the previous known location. 
'''

def ExtractBHLocations(source,BH1PrevLoc=[0.0,0.0,0.0],BH2PrevLoc=[0.0,0.0,0.0],varName='U_CHI',tol=1e-3):
    polyData=source.GetOutput()
    pointData=polyData.GetPointData()
    data=pointData.GetArray(varName)
    data_range=data.GetRange()
    minValue=data_range[0]
    numPoints=polyData.GetNumberOfPoints()
    bhLoc=[[0,0,0],[0,0,0]]

    BH1PrevLoc=np.array(BH1PrevLoc)
    BH2PrevLoc=np.array(BH2PrevLoc)
    BH1CurrLoc=np.zeros((1,3))
    BH2CurrLoc=np.zeros((1,3))

    BH1L2=float('inf')
    BH2L2=float('inf')

    for p in range(numPoints):
        value=data.GetComponent(p,0)
        if abs(value-minValue)<tol:
            # valid candidate
            validMin=[0,0,0]
            polyData.GetPoint(p,validMin)
            validMin=np.array(validMin)
            distBH1=np.linalg.norm(validMin-BH1PrevLoc)
            distBH2=np.linalg.norm(validMin-BH2PrevLoc)
            
            if(distBH1<distBH2):
                if(np.linalg.norm(validMin-BH1PrevLoc)<BH1L2):
                    BH1L2=np.linalg.norm(validMin-BH1PrevLoc)
                    BH1CurrLoc=validMin
            elif(distBH1>distBH2):
                if(np.linalg.norm(validMin-BH2PrevLoc)<BH2L2):
                    BH2L2=np.linalg.norm(validMin-BH2PrevLoc)
                    BH2CurrLoc=validMin
            else:
                if(np.linalg.norm(validMin-BH1PrevLoc)<BH1L2):
                    BH1L2=np.linalg.norm(validMin-BH1PrevLoc)
                    BH1CurrLoc=validMin

                if(np.linalg.norm(validMin-BH2PrevLoc)<BH2L2):
                    BH2L2=np.linalg.norm(validMin-BH2PrevLoc)
                    BH2CurrLoc=validMin


            

            
    bhLoc=[BH1CurrLoc.tolist(),BH2CurrLoc.tolist()]
    return bhLoc






