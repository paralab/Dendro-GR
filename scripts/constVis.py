
# coding: utf-8

# In[3]:


# Script to plot Psi4 data
import argparse
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

parser = argparse.ArgumentParser(description='generate seqence of images')
parser.add_argument('-b','--begin', help='begin step', required=True)
parser.add_argument('-e','--end', help='end step', required=True)
parser.add_argument('-f','--freq', help='io step frequency', required=True)
parser.add_argument('-pvtu','--pvtu_prefix', help='pvtu prefix', required=True)
parser.add_argument('-img','--img_prefix', help='image prefix ', required=True)
args = vars(parser.parse_args())


# In[11]:


def saveImg(file_prefix,img_prefix,step,imgIndex):

    filename=[file_prefix+'_'+str(step)+'.pvtu']
    print "reading: %s" %(filename)
    # create a new 'XML Partitioned Unstructured Grid Reader'
    bssn_gr_const_0pvtu = XMLPartitionedUnstructuredGridReader(FileName=filename)
    bssn_gr_const_0pvtu.CellArrayStatus = ['mpi_rank', 'cell_level']
    bssn_gr_const_0pvtu.PointArrayStatus = ['C_HAM', 'C_MOM0', 'C_MOM1', 'C_MOM2', 'C_PSI4_REAL', 'C_PSI4_IMG']

    # Properties modified on bssn_gr_0pvtu
    bssn_gr_const_0pvtu.CellArrayStatus = ['cell_level']
    bssn_gr_const_0pvtu.PointArrayStatus = ['C_PSI4_REAL']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1176, 808]
    #renderView1.ViewSize = [4096 , 2160]
    #renderView1.ViewSize = [1600,900]
    #renderView1.ViewSize = [800,600]
    
    #print 'reaach 1'
    # show data in view
    bssn_gr_const_0pvtuDisplay = Show(bssn_gr_const_0pvtu, renderView1)
    # trace defaults for the display properties.
    bssn_gr_const_0pvtuDisplay.ColorArrayName = [None, '']
    bssn_gr_const_0pvtuDisplay.ScalarOpacityUnitDistance = 106.21591469472322

    # reset view to fit data
    renderView1.ResetCamera()
    #renderView1.Zoom(2.0)
    renderView1.CameraPosition = [2048.0, 2048.0, 2350.8216648370403]
    #renderView1.CameraPosition = [2048.0, 2048.0, 3150.0]
    renderView1.CameraFocalPoint = [2048.0, 2048.0, 2048.0]
    renderView1.CameraParallelScale = 3500


    # create a new 'Slice'
    slice1 = Slice(Input=bssn_gr_const_0pvtu)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [2048.0, 2048.0, 2048.0]
    #slice1.SliceType.Origin = [1024.0, 1024.0, 1024.0]

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1)
    # trace defaults for the display properties.
    slice1Display.ColorArrayName = [None, '']

    #print 'reaach 2'
    # hide data in view
    Hide(bssn_gr_const_0pvtu, renderView1)

    # set scalar coloring
    ColorBy(slice1Display, ('POINTS', 'C_PSI4_REAL'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'UCHI'
    uCHILUT = GetColorTransferFunction('CPSI4REAL')

    # get opacity transfer function/opacity map for 'UCHI'
    uCHIPWF = GetOpacityTransferFunction('CPSI4REAL')

    # Properties modified on uCHILUT
    uCHILUT.ColorSpace = 'HSV'

    # Properties modified on uCHILUT
    #uCHILUT.LockScalarRange = 1

    # Properties modified on uCHILUT
    uCHILUT.NumberOfTableValues = 1024
    annotateTimeFilter1 = AnnotateTimeFilter(Input=slice1Display)
    # Properties modified on annotateTimeFilter1
    annotateTimeFilter1.Format = 'Time: %f s' %(0.1*(400.0/(1<<12))*float(step))
    annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

    #print 'reaach 2.5'
    # export view
    imgName=img_prefix+'_slice_'+str(imgIndex).zfill(6)+'.png'
    viewLayout1 = GetLayout()
    SaveScreenshot(imgName, layout=viewLayout1, magnification=1, quality=100)
    #ExportView(imgName, view=renderView1)
    # Properties modified on uCHILUT

    #print 'reaach 3'
    slice1Display.SetScalarBarVisibility(renderView1, False)
    # set scalar coloring
    ColorBy(slice1Display, ('CELLS', 'cell_level'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'celllevel'
    celllevelLUT = GetColorTransferFunction('celllevel')

    # get opacity transfer function/opacity map for 'celllevel'
    celllevelPWF = GetOpacityTransferFunction('celllevel')

    #print 'reaach 3.5'
    # export view
    #ExportView('/home/milinda/Desktop/gr_sp1/slice_level.svg', view=renderView1)
    imgName=img_prefix+'_slice_level_'+str(imgIndex).zfill(6)+'.png'
    #ExportView(imgName, view=renderView1)
    viewLayout1 = GetLayout()
    SaveScreenshot(imgName, layout=viewLayout1, magnification=1, quality=100)

    slice1Display.SetScalarBarVisibility(renderView1, False)



    # hide data in view
    #Hide(warpByScalar1, renderView1)

    # set active source
    #SetActiveSource(slice1)

    Delete(annotateTimeFilter1Display)
    del annotateTimeFilter1Display

    Delete(annotateTimeFilter1)
    del annotateTimeFilter1

    # set active source
    SetActiveSource(bssn_gr_const_0pvtu)

    # destroy slice1
    Delete(slice1)
    del slice1

    # destroy bssn_gr_0pvtu
    Delete(bssn_gr_const_0pvtu)
    del bssn_gr_const_0pvtu

    # destroy renderView1
    Delete(renderView1)
    del renderView1

# In[12]:


step_begin=int(args['begin'])
step_end=int(args['end'])
step_freq=int(args['freq'])
file_prefix=args['pvtu_prefix']
img_prefix=args['img_prefix']


#step_begin=0
#step_end=2000
#step_freq=50
#file_prefix='/home/milinda/Desktop/gr_sp1'

n_procs = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
print "number of ranks: ", n_procs
imgCount=step_begin/step_freq;
for step in range(step_begin,step_end,step_freq):
    # create a new 'XML Partitioned Unstructured Grid Reader'
    saveImg(file_prefix,img_prefix,step,imgCount)
    imgCount=imgCount+1;

#
