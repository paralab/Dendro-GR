# author: Milinda Fernando
# simple paraview python script to bh binary visualization problem ().
# School of Computing, University of Utah
# 12/12/2017
import argparse
import resource
import gc
#from pympler.tracker import SummaryTracker
# import the simple module from the paraview
from paraview.simple import *
#import vtk
import paraview.vtk.numpy_interface.dataset_adapter as dsa
import paraview.vtk.numpy_interface.algorithms as algs
import json as js
from pprint import pprint
#from paraview.simple import vtk

# disable automatic camera reset on 'Show'

paraview.simple._DisableFirstRenderCameraReset()


def filter_warp_by_scalar(inp, scale, varName):
    warpByScalar = WarpByScalar(Input=inp)
    warpByScalar.Scalars = ['POINTS', varName]
    warpByScalar.ScaleFactor=scale
    #WbySDisplay.ColorArrayName = [None, '']
    #WbySDisplay.SetRepresentationType(params["PARAVIEW_DISPLAY_MODE"])
    #ColorBy(WbySDisplay, ('POINTS', varName))
    return warpByScalar

def filter_threshold(inp,varType,varName,rng):
    threshold = Threshold(Input=inp)
    threshold.Scalars = [varType, varName]
    threshold.ThresholdRange = [rng[0],rng[1]]
    return threshold


def filter_slice(inp,sOrigin,sNormal,sliceType='Plane'):
     # create a new 'Slice'
    slice1 = Slice(Input=inp)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = sOrigin
    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = sNormal
    slice1.Triangulatetheslice = 0

    return slice1
'''
computes and add r \times \psi4 for the input data
'''
def add_rxpsi4(inp,params,vname):

    cal1 = Calculator(Input=inp)
    # Properties modified on calculator1
    cal1.ResultArrayName = 'posVec'
    cal1.Function = '(-200 + coordsX*(400.0/2^14))*iHat + (-200 + coordsY*(400.0/2^14))*jHat + (-200 + coordsZ*(400.0/2^14))*kHat'

    cal2 = Calculator(Input=cal1)
    cal2.ResultArrayName = vname
    cal2.Function = 'mag(posVec)*sqrt(C_PSI4_REAL^2 + C_PSI4_IMG^2)'

    return cal2
    

def saveImg(params, step, imgIndex,cam_coords):

    maxDepth = params["BSSN_MAXDEPTH"]
    CHI_FLOOR = params["CHI_FLOOR"]

    _evolveVars = ["U_ALPHA", "U_CHI", "U_K", "U_GT0", "U_GT1", "U_GT2", "U_BETA0", "U_BETA1", "U_BETA2", "U_B0", "U_B1", "U_B2", "U_SYMGT0",
                   "U_SYMGT1", "U_SYMGT2", "U_SYMGT3", "U_SYMGT4", "U_SYMGT5", "U_SYMAT0", "U_SYMAT1", "U_SYMAT2", "U_SYMAT3", "U_SYMAT4", "U_SYMAT5"]
    _constVars = ["C_HAM", "C_MOM0", "C_MOM1",
                  "C_MOM2", "C_PSI4_REAL", "C_PSI4_IMG","C_R_PSI4"]

    file_prefix = params["BSSN_VTU_FILE_PREFIX"]
    img_prefix = params["BSSN_IMG_FILE_PREFIX"]

    s_x = (1 << (maxDepth-1))
    s_y = (1 << (maxDepth-1))
    s_z = (1 << (maxDepth-1))

    scale_factor = (1 << (maxDepth))
    #tracker = SummaryTracker()
    # mem_begin=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    filename = [str(file_prefix+"_"+str(step)+".pvtu")]
    print ("reading: %s" % (filename))
    # create a new 'XML Partitioned Unstructured Grid Reader'

    bssn_gr_0pvtu = XMLPartitionedUnstructuredGridReader(FileName=filename)

    numEVolVars = params["BSSN_NUM_EVOL_VARS_VTU_OUTPUT"]
    numConstVars = params["BSSN_NUM_CONST_VARS_VTU_OUTPUT"]
    numTotal = numEVolVars+numConstVars

    varNames = []

    for index in params["BSSN_VTU_OUTPUT_EVOL_INDICES"]:
        varNames.append(_evolveVars[index])

    for index in params["BSSN_VTU_OUTPUT_CONST_INDICES"]:
        varNames.append(_constVars[index])

    # print varNames

    # Properties modified on bssn_gr_0pvtu
    bssn_gr_0pvtu.CellArrayStatus = ['cell_level']
    bssn_gr_0pvtu.PointArrayStatus = varNames

    bssn_gr_0pvtu=add_rxpsi4(bssn_gr_0pvtu,params,'C_R_PSI4')

    UpdatePipeline()
    print ("reading ended")
    

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [params["VTK_RENDER_VIEW_SIZE"]
                            [0], params["VTK_RENDER_VIEW_SIZE"][1]]

    # print 'reaach 1'
    bssn_gr_0pvtuDisplay = Show(bssn_gr_0pvtu, renderView1)
    # trace defaults for the display properties.
    bssn_gr_0pvtuDisplay.ColorArrayName = [None, '']
    bssn_gr_0pvtuDisplay.ScalarOpacityUnitDistance = 106.21591469472322

    # reset view to fit data
    renderView1.ResetCamera()
    if len(cam_coords)==0:
        renderView1.CameraPosition = [s_x, s_y, 0.7*s_z]
    else:
        renderView1.CameraPosition = [cam_coords[0], cam_coords[1], cam_coords[2]]
    renderView1.CameraFocalPoint = [s_x, s_y, s_z]
    renderView1.CameraParallelScale = 3500

    # UpdatePipeline()
    renderView1.Update()

    slice1=filter_slice(bssn_gr_0pvtu,sOrigin=[s_x,s_y,s_z],sNormal=[0,0,1])
    Hide(bssn_gr_0pvtu,renderView1)
    print ("slice created")
    UpdatePipeline()
    print ("pipeline updated")

    for varName in varNames:
        # show data in view
        wbs=filter_warp_by_scalar(slice1,-300,varName)
        wbsDisplay=Show(wbs,renderView1)
        
        Hide(bssn_gr_0pvtu, renderView1)
        Hide(slice1, renderView1)
        UpdatePipeline()

        activeDisplay=wbsDisplay
        
        # trace defaults for the display properties.
        #slice1Display.ColorArrayName = [None, '']
        #slice1Display.SetRepresentationType(params["PARAVIEW_DISPLAY_MODE"])

        activeDisplay.ColorArrayName = [None, '']
        activeDisplay.SetRepresentationType(params["PARAVIEW_DISPLAY_MODE"])

       

        renderView1.OrientationAxesVisibility = 0
        ColorBy(activeDisplay, ('POINTS', varName))

        # get color transfer function/color map for 'UCHI'
        uCHILUT = GetColorTransferFunction(varName)

        # get opacity transfer function/opacity map for 'UCHI'
        uCHIPWF = GetOpacityTransferFunction(varName)

        # Properties modified on uCHILUT
        uCHILUT.ColorSpace = 'Diverging'
        uCHILUT.NumberOfTableValues = 1024
        u_CHILUTColorBar = GetScalarBar(uCHILUT, renderView1)
        u_CHILUTColorBar.Title = varName
        u_CHILUTColorBar.ComponentTitle = ''
        u_CHILUTColorBar.TitleJustification = 'Centered'

        #uCHILUT.RescaleTransferFunction(CHI_FLOOR, 1.0)
        #uCHIPWF.RescaleTransferFunction(CHI_FLOOR, 1.0)
        activeDisplay.RescaleTransferFunctionToDataRange(False, True)

        if(params["SET_COLOR_BAR_VISIBILITY"] == 1):
            activeDisplay.SetScalarBarVisibility(renderView1, True)
        else:
            activeDisplay.SetScalarBarVisibility(renderView1, True)

        UpdatePipeline()
        renderView1.Update()

        bh_range_max = params["BSSN_GRID_MAX_X"]
        bh_range_min = params["BSSN_GRID_MIN_X"]

        if(params["ANNOTATE_TIME"] == 1):
            annotateTimeFilter1 = AnnotateTimeFilter(Input=activeDisplay)
            # Properties modified on annotateTimeFilter1
            annotateTimeFilter1.Format = 'Time: %f s' % (
                params["BSSN_CFL_FACTOR"]*(float(bh_range_max-bh_range_min)/(1 << maxDepth))*float(step))
            annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)
            UpdatePipeline()

        # export view
        imgName = str(img_prefix+"_slice_"+str(varName) +
                      "_"+str(imgIndex).zfill(6)+".png")
        #viewLayout1 = GetLayout()
        SaveScreenshot(imgName, magnification=1, quality=100, view=renderView1)
        #ExportView(imgName, view=renderView1)
        print ("image written for variable "+str(varName))

        if(params["SET_COLOR_BAR_VISIBILITY"] == 1):
            activeDisplay.SetScalarBarVisibility(renderView1, False)

        SetActiveSource(wbs)
        Delete(wbs)
        #del wbs

        #UpdatePipeline()
        renderView1.Update()

    # set active source
    # destroy slice1
    SetActiveSource(slice1)
    Delete(slice1)
    del slice1

    # destroy bssn_gr_0pvtu
    SetActiveSource(bssn_gr_0pvtu)
    Delete(bssn_gr_0pvtu)
    del bssn_gr_0pvtu

    # destroy renderView1
    Delete(renderView1)
    del renderView1


def writeZSlice(params):

    maxDepth = params["BSSN_MAXDEPTH"]
    CHI_FLOOR = params["CHI_FLOOR"]

    _evolveVars = ["U_ALPHA", "U_CHI", "U_K", "U_GT0", "U_GT1", "U_GT2", "U_BETA0", "U_BETA1", "U_BETA2", "U_B0", "U_B1", "U_B2", "U_SYMGT0",
                   "U_SYMGT1", "U_SYMGT2", "U_SYMGT3", "U_SYMGT4", "U_SYMGT5", "U_SYMAT0", "U_SYMAT1", "U_SYMAT2", "U_SYMAT3", "U_SYMAT4", "U_SYMAT5"]
    _constVars = ["C_HAM", "C_MOM0", "C_MOM1",
                  "C_MOM2", "C_PSI4_REAL", "C_PSI4_IMG","C_R_PSI4"]

    file_prefix = params["BSSN_VTU_FILE_PREFIX"]
    img_prefix = params["BSSN_IMG_FILE_PREFIX"]

    s_x = (1 << (maxDepth-1))
    s_y = (1 << (maxDepth-1))
    s_z = (1 << (maxDepth-1))

    scale_factor = (1 << (maxDepth))
    #tracker = SummaryTracker()
    # mem_begin=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    filename = list()
    print("reading: %s" % (filename))
    # create a new 'XML Partitioned Unstructured Grid Reader'

    for step in range(params["BSSN_TIMESTEP_BEGIN"], params["BSSN_TIMESTEP_END"], params["BSSN_IO_OUTPUT_FREQ"]):
        filename.append(str(file_prefix+"_"+str(step)+".pvtu"))

    bssn_gr_0pvtu = XMLPartitionedUnstructuredGridReader(FileName=filename)

    #renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    #renderView1.ViewSize = [params["VTK_RENDER_VIEW_SIZE"]
    #                        [0], params["VTK_RENDER_VIEW_SIZE"][1]]

    # print 'reaach 1'
    #bssn_gr_0pvtuDisplay = Show(bssn_gr_0pvtu, renderView1)

    # get animation scene
    #animationScene1 = GetAnimationScene()

    # update animation scene based on data timesteps
    #animationScene1.UpdateAnimationUsingDataTimeSteps()

    numEVolVars = params["BSSN_NUM_EVOL_VARS_VTU_OUTPUT"]
    numConstVars = params["BSSN_NUM_CONST_VARS_VTU_OUTPUT"]
    numTotal = numEVolVars+numConstVars

    varNames = []

    for index in params["BSSN_VTU_OUTPUT_EVOL_INDICES"]:
        varNames.append(_evolveVars[index])

    for index in params["BSSN_VTU_OUTPUT_CONST_INDICES"]:
        varNames.append(_constVars[index])

    # print varNames

    # Properties modified on bssn_gr_0pvtu
    bssn_gr_0pvtu.CellArrayStatus = ['cell_level']
    bssn_gr_0pvtu.PointArrayStatus = varNames
    #bssn_gr_0pvtu.CellArrayStatus = ['mpi_rank', 'cell_level']
    #bssn_gr_0pvtu.PointArrayStatus = ['U(x,t)']


    #bssn_gr_0pvtu=add_rxpsi4(bssn_gr_0pvtu,params,'C_R_PSI4')

    UpdatePipeline()
    print ("reading ended")

    slice1=filter_slice(bssn_gr_0pvtu,sOrigin=[s_x,s_y,s_z],sNormal=[0,0,1])
    print ("slice created")

    #slice1Display = Show(slice1, renderView1)

    # save data
    #SaveData(file_prefix+"_test.pvd", proxy=slice1,Writetimestepsasfileseries=1, Firsttimestep=0, Lasttimestep=params["BSSN_TIMESTEP_END"]-params["BSSN_TIMESTEP_BEGIN"], Timestepstride=1, Filenamesuffix='_%d', DataMode='Appended', CompressorType='Zlib')
    fileslice = file_prefix+'_slice.pvd'
    SaveData(fileslice, proxy=slice1, CompressorType='LZ4', WriteTimeSteps=1)
    '''
    SaveData(fileslice, proxy=slice1, Writetimestepsasfileseries=1,
    Firsttimestep=int(params["BSSN_TIMESTEP_BEGIN"]),
    Lasttimestep=int(params["BSSN_TIMESTEP_END"]),
    Timestepstride=int(params["BSSN_IO_OUTPUT_FREQ"]),
    Filenamesuffix='_%d',
    DataMode='Appended',
    EncodeAppendedData=0,
    CompressorType='None',
    CompressionLevel='6') '''
    #for step in range(params["BSSN_TIMESTEP_BEGIN"], params["BSSN_TIMESTEP_END"], params["BSSN_IO_OUTPUT_FREQ"]):
        
    #    filename_slice = str(file_prefix+"_"+str(step)+".pvd")
        
     
        #SaveData( filename_slice[0] , proxy=slice1, CompressorType='ZLib')
    




def main():
    # print command line arguments
    if(len(sys.argv) == 0):
        print ("Error: Visualization parameter file not specified")
        sys.exit(0)

    with open(sys.argv[1]) as f:
        params = js.load(f)

    #print("parameters read")
    # pprint(params)
    
    maxDepth = params["BSSN_MAXDEPTH"]
    s_x = (1 << (maxDepth-1))
    s_y = (1 << (maxDepth-1))
    s_z = (1 << (maxDepth-1))

    zfac=0.2*((params["BSSN_TIMESTEP_END"]-params["BSSN_TIMESTEP_BEGIN"])/(float)(params["BSSN_IO_OUTPUT_FREQ"]))
    zfac=(0.8*s_z + 2*s_z)/zfac
    
    yfac=0.2*((params["BSSN_TIMESTEP_END"]-params["BSSN_TIMESTEP_BEGIN"])/(float)(params["BSSN_IO_OUTPUT_FREQ"]))
    yfac=(2*s_y)/yfac
    #print(zfac)

    params["BSSN_TIMESTEP_BEGIN"] = int(sys.argv[2])
    params["BSSN_TIMESTEP_END"] = int(sys.argv[3])

    n_procs = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
    print ("number of ranks: ", n_procs)
    imgCount = params["BSSN_TIMESTEP_BEGIN"]/params["BSSN_IO_OUTPUT_FREQ"]
    #writeZSlice(params)
    for step in range(params["BSSN_TIMESTEP_BEGIN"], params["BSSN_TIMESTEP_END"], params["BSSN_IO_OUTPUT_FREQ"]):
        saveImg(params, step, imgCount,[s_x,min(-s_y + imgCount*yfac,s_y), min(-2*s_z+(imgCount*zfac),0.8*s_z)])
        #writeZSlice(params,step,imgCount)
        imgCount = imgCount+1

#


if __name__ == "__main__":
    main()
