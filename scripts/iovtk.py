##
# @author: Milinda Fernando
# School of Computing, University of Utah. 
# @brief: Contains code for reading the vtk data files. 
# @date: 05/05/2019

import os
import vtk as vtk
import numpy as np
from mpi4py import MPI
from vtk.util import numpy_support


'''
@brief Reads the XML unstructed grid and retuns, vtkXMLUnstructuredGridReader (sequential)
'''
def ReadVTUFile(fname):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    return reader    

def WriteVTUFile(fname, source):
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(fname)
    #ug = vtk.vtkUnstructuredGrid()
    #ug.setPoints(source.GetPoints())
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection( source.GetOutputPort() )
    
    writer.SetInputData(mapper.GetInput())
    writer.Write()
    
'''
@brief Reads the XML partioned unstructed grid and retuns, vtkXMLPUnstructuredGridReader
'''
def ReadPVTUFile(fname):
    reader = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()

    return reader

'''
@brief : convert vtk to numpy data. 
'''
def vtk_to_numpy(vtkReader,varName):
    vtk_data   = vtkReader.GetOutput().GetPointData()
    raw_data =  vtk_data.GetArray(varName)
    #print(raw_data.GetNumberOfTuples())
    #print(raw_data.GetNumberOfComponents())
    #print(raw_data.GetCol(0))
    #print(len(numpy_data[0]))
    #print(numpy_data)
    #print(numpy_data[0][2])

    numpy_data  = numpy_support.vtk_to_numpy(raw_data)
    v=list()
    if(raw_data.GetNumberOfComponents() > 1):
        for i in range(0,raw_data.GetNumberOfComponents()):
            v.append(numpy_data[:,i])
        return v
    else:
        return numpy_data

    
    


'''
@brief: convert vtk to numpy mat, based on the number of timesteps to read. 
'''
def vtk_to_numpy_mat(fprefix,tbegin,tend,varName):
    
    fname = str(fprefix) + str(tbegin).zfill(4) + ".pvtu"
    reader = ReadPVTUFile(fname)
    raw_data=reader.GetOutput().GetPointData().GetArray(varName)

    numC = raw_data.GetNumberOfComponents()
    n = raw_data.GetNumberOfTuples()
    T = (tend-tbegin)
    #print(numC)
    if(numC > 1):
        numpy_data=[]
        for i in range(0,numC):
            comp_i = np.zeros((T,n))
            numpy_data.append(comp_i)
        
        for step in range(tbegin,tend):
            fname = str(fprefix) + str(step).zfill(4) + ".pvtu"
            reader = ReadPVTUFile(fname)
            np_row = vtk_to_numpy(reader,varName)
            for i in range(0,numC):
                numpy_data[i][(step-tbegin)] = np_row[i]

    else:
        numpy_data = np.zeros((T,n))
        for step in range(tbegin,tend):
            fname = str(fprefix) + str(step).zfill(4) + ".pvtu"
            reader = ReadPVTUFile(fname)
            np_row = vtk_to_numpy(reader,varName)
            numpy_data[(step-tbegin)] = np_row

    return numpy_data
