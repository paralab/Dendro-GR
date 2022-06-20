
# coding: utf-8

# In[3]:


# author: Milinda Fernando
# simple paraview python script to bh binary visualization problem.
# School of Computing, University of Utah
# 12/12/2017
import argparse
import resource
import gc
#from pympler.tracker import SummaryTracker
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

    maxDepth=12;
    bh_range_min=-200
    bh_range_max=200
    s_x=(1<<(maxDepth-1))
    s_y=(1<<(maxDepth-1))
    s_z=(1<<(maxDepth-1))
    scale_factor=2048.0;
    #tracker = SummaryTracker()
    #mem_begin=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    filename=[file_prefix+'_'+str(step)+'.pvtu']
    print "reading: %s" %(filename)
    # create a new 'XML Partitioned Unstructured Grid Reader'
    bssn_gr_0pvtu = XMLPartitionedUnstructuredGridReader(FileName=filename)
    bssn_gr_0pvtu.CellArrayStatus = ['mpi_rank', 'cell_level']
    bssn_gr_0pvtu.PointArrayStatus = ['U_ALPHA', 'U_CHI', 'U_K', 'U_GT0', 'U_GT1', 'U_GT2', 'U_BETA0', 'U_BETA1', 'U_BETA2', 'U_B0', 'U_B1', 'U_B2', 'U_SYMGT0', 'U_SYMGT1', 'U_SYMGT2', 'U_SYMGT3', 'U_SYMGT4', 'U_SYMGT5', 'U_SYMAT0', 'U_SYMAT1', 'U_SYMAT2', 'U_SYMAT3', 'U_SYMAT4', 'U_SYMAT5']

    # Properties modified on bssn_gr_0pvtu
    bssn_gr_0pvtu.CellArrayStatus = ['cell_level']
    bssn_gr_0pvtu.PointArrayStatus = ['U_ALPHA', 'U_CHI', 'U_K']

    UpdatePipeline()
    print "reading ended"

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    #renderView1.ViewSize = [1176, 808]
    #renderView1.ViewSize = [4096 , 2160]
    renderView1.ViewSize = [1600,900]

    #print 'reaach 1'
    # show data in view
    bssn_gr_0pvtuDisplay = Show(bssn_gr_0pvtu, renderView1)
    # trace defaults for the display properties.
    bssn_gr_0pvtuDisplay.ColorArrayName = [None, '']
    bssn_gr_0pvtuDisplay.ScalarOpacityUnitDistance = 106.21591469472322

    # reset view to fit data
    renderView1.ResetCamera()
    #renderView1.Zoom(2.0)
    #renderView1.CameraPosition = [2048.0, 2048.0, 2350.8216648370403]
    renderView1.CameraPosition = [s_x, s_y, 3150.0*s_z/scale_factor]
    renderView1.CameraFocalPoint = [s_x, s_y, s_z]
    renderView1.CameraParallelScale = 3500

    programmableFilter1 = ProgrammableFilter(Input=bssn_gr_0pvtu)
    programmableFilter1.Script = ''
    programmableFilter1.RequestInformationScript = ''
    programmableFilter1.RequestUpdateExtentScript = ''
    programmableFilter1.PythonPath = ''

    # Properties modified on programmableFilter1
    with open('trajectory_filter.py', 'r') as filter_trajectory:
        data=filter_trajectory.read()
        filter_trajectory.close()


    programmableFilter1.Script =data #'input=inputs[0]\n u_chi=input.PointData[\'U_CHI\']\nnumPoints=input.GetNumberOfPoints()\nu_chi_min=min(u_chi)#u_chi[0] \nfor i in range(numPoints):\n    uchi=u_chi[i]\n    if abs(uchi-u_chi_min)< 1e-5 :\n        minpos=input.GetPoint(i)\n        #u_chi_min=uchi\n        print "uzmin: ",uchi,min(u_chi)\n        print "minpos: ",minpos\n        break;\n        #output.append(u_chi_min, \'u_chi_min\')\n        #output.append(u_chi_min_pos, \'u_chi_min_pos\')```    '
    programmableFilter1.RequestInformationScript = ''
    programmableFilter1.RequestUpdateExtentScript = ''
    programmableFilter1.PythonPath = ''
    programmableFilter1Display = Show(programmableFilter1, renderView1)

    SetActiveSource(programmableFilter1)
    Delete(programmableFilter1)
    del programmableFilter1

    #UpdatePipeline()
    renderView1.Update()



    # create a new 'Slice'
    slice1 = Slice(Input=bssn_gr_0pvtu)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [s_x, s_y, s_z]

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
    slice1.Triangulatetheslice = 0
    
    print "slice created"
    UpdatePipeline()
    print "pipeline updated"
    # show data in view
    slice1Display = Show(slice1, renderView1)
    # trace defaults for the display properties.
    slice1Display.ColorArrayName = [None, '']

    #print 'reaach 95'
    # hide data in view
    Hide(bssn_gr_0pvtu, renderView1)
    UpdatePipeline()
    renderView1.OrientationAxesVisibility = 0
    ColorBy(slice1Display, ('POINTS','U_CHI'))
    slice1Display.SetScalarBarVisibility(renderView1, True)
    #print 'reaach 99'
    # set scalar coloring
    # get color transfer function/color map for 'U_CHI'
    u_CHILUT = GetColorTransferFunction('U_CHI')
    u_CHILUT.LockDataRange = 0
    u_CHILUT.InterpretValuesAsCategories = 0
    u_CHILUT.ShowCategoricalColorsinDataRangeOnly = 0
    u_CHILUT.RescaleOnVisibilityChange = 0
    u_CHILUT.EnableOpacityMapping = 0
    u_CHILUT.RGBPoints = [6.286501118071329e-08, 0.231373, 0.298039, 0.752941, 0.4964692923160348, 0.865003, 0.865003, 0.865003, 0.9929385217670584, 0.705882, 0.0156863, 0.14902]
    u_CHILUT.UseLogScale = 0
    u_CHILUT.ColorSpace = 'Diverging'
    u_CHILUT.UseBelowRangeColor = 0
    u_CHILUT.BelowRangeColor = [0.0, 0.0, 0.0]
    u_CHILUT.UseAboveRangeColor = 0
    u_CHILUT.AboveRangeColor = [1.0, 1.0, 1.0]
    u_CHILUT.NanColor = [1.0, 1.0, 0.0]
    u_CHILUT.Discretize = 1
    u_CHILUT.NumberOfTableValues = 1024
    u_CHILUT.ScalarRangeInitialized = 1.0
    u_CHILUT.HSVWrap = 0
    u_CHILUT.VectorComponent = 0
    u_CHILUT.VectorMode = 'Magnitude'
    u_CHILUT.AllowDuplicateScalars = 1
    u_CHILUT.Annotations = []
    u_CHILUT.ActiveAnnotatedValues = []
    u_CHILUT.IndexedColors = []

    # toggle 3D widget visibility (only when running from the GUI)
    #Hide3DWidgets(proxy=slice1.SliceType)

    # get color legend/bar for u_CHILUT in view renderView1
    u_CHILUTColorBar = GetScalarBar(u_CHILUT, renderView1)
    u_CHILUTColorBar.AutoOrient = 1
    u_CHILUTColorBar.Orientation = 'Vertical'
    u_CHILUTColorBar.WindowLocation = 'LowerRightCorner'
    u_CHILUTColorBar.Position = [0.89, 0.02]
    u_CHILUTColorBar.Title = "U_CHI"
    u_CHILUTColorBar.ComponentTitle = ''
    u_CHILUTColorBar.TitleJustification = 'Centered'
    u_CHILUTColorBar.TitleColor = [1.0, 1.0, 1.0]
    u_CHILUTColorBar.TitleOpacity = 1.0
    u_CHILUTColorBar.TitleFontFamily = 'Arial'
    u_CHILUTColorBar.TitleBold = 0
    u_CHILUTColorBar.TitleItalic = 0
    u_CHILUTColorBar.TitleShadow = 0
    u_CHILUTColorBar.TitleFontSize = 16
    u_CHILUTColorBar.LabelColor = [1.0, 1.0, 1.0]
    u_CHILUTColorBar.LabelOpacity = 1.0
    u_CHILUTColorBar.LabelFontFamily = 'Arial'
    u_CHILUTColorBar.LabelBold = 0
    u_CHILUTColorBar.LabelItalic = 0
    u_CHILUTColorBar.LabelShadow = 0
    u_CHILUTColorBar.LabelFontSize = 16
    u_CHILUTColorBar.AutomaticLabelFormat = 0
    u_CHILUTColorBar.LabelFormat = '%-#2.2g'
    u_CHILUTColorBar.DrawTickMarks = 1
    u_CHILUTColorBar.DrawTickLabels = 1
    u_CHILUTColorBar.UseCustomLabels = 0
    u_CHILUTColorBar.CustomLabels = []
    u_CHILUTColorBar.AddRangeLabels = 1
    u_CHILUTColorBar.RangeLabelFormat = '%-#6.1e'
    u_CHILUTColorBar.DrawAnnotations = 1
    u_CHILUTColorBar.AddRangeAnnotations = 0
    u_CHILUTColorBar.AutomaticAnnotations = 0
    u_CHILUTColorBar.DrawNanAnnotation = 0
    u_CHILUTColorBar.NanAnnotation = 'NaN'
    u_CHILUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
    u_CHILUTColorBar.ScalarBarThickness = 16
    u_CHILUTColorBar.ScalarBarLength = 0.33

    # Rescale transfer function

    # get opacity transfer function/opacity map for 'U_CHI'
    u_CHIPWF = GetOpacityTransferFunction('U_CHI')
    u_CHIPWF.Points = [6.286501118071329e-08, 0.0, 0.5, 0.0, 0.9929385217670584, 1.0, 0.5, 0.0]
    u_CHIPWF.AllowDuplicateScalars = 1
    u_CHIPWF.ScalarRangeInitialized = 1

    # Rescale transfer function
    u_CHILUT.RescaleTransferFunction(0.0001, 1.0)
    u_CHIPWF.RescaleTransferFunction(0.0001, 1.0)
    UpdatePipeline()
    renderView1.Update()
    #annotateTimeFilter1 = AnnotateTimeFilter(Input=slice1Display)
    # Properties modified on annotateTimeFilter1
    #annotateTimeFilter1.Format = 'Time: %f s' %(0.1*(float(bh_range_max-bh_range_min)/(1<<maxDepth))*float(step))
    #annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

    #print 'reaach 129'
    #print 'reaach 2.5'
    # export view
    imgName=img_prefix+'_slice_'+str(imgIndex).zfill(6)+'.png'
    #viewLayout1 = GetLayout()
    #print 'reaach 134'
    SaveScreenshot(imgName,magnification=1, quality=100,view=renderView1)
    #print "image1 written"
    #Delete(viewLayout1)
    #del viewLayout1
    #ExportView(imgName, view=renderView1)
    # Properties modified on uCHILUT

    #print 'reaach 3'
    slice1Display.SetScalarBarVisibility(renderView1, False)
    # set scalar coloring
    ColorBy(slice1Display, ('CELLS', 'cell_level'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True,False)



    # get color transfer function/color map for 'celllevel'
    celllevelLUT = GetColorTransferFunction('cell_level')

    # get opacity transfer function/opacity map for 'celllevel'
    celllevelPWF = GetOpacityTransferFunction('cell_level')
    celllevelLUTColorBar = GetScalarBar(celllevelLUT, renderView1)
    celllevelLUTColorBar.Title="level"
    celllevelLUTColorBar.ComponentTitle=''
    celllevelLUTColorBar.AutomaticLabelFormat = 0
    celllevelLUTColorBar.LabelFormat = '%-#.2f'
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    UpdatePipeline()
    renderView1.Update()
    #print 'reaach 3.5'
    # export view
    #ExportView('/home/milinda/Desktop/gr_sp1/slice_level.svg', view=renderView1)
    imgName=img_prefix+'_slice_level_'+str(imgIndex).zfill(6)+'.png'
    #ExportView(imgName, view=renderView1)
    #viewLayout1 = GetLayout()
    SaveScreenshot(imgName,magnification=1, quality=100,view=renderView1)
    #Delete(viewLayout1)
    #del viewLayout1
    print "image2 written"

    slice1Display.SetScalarBarVisibility(renderView1, False)
    # create a new 'Warp By Scalar'
    warpByScalar1 = WarpByScalar(Input=slice1)
    #warpByScalar1.Scalars = ['POINTS', 'U_ALPHA']

    # Properties modified on warpByScalar1
    warpByScalar1.Scalars = ['POINTS', 'U_CHI']
    warpByScalar1.ScaleFactor = 500.0

    # show data in view
    warpByScalar1Display = Show(warpByScalar1, renderView1)
    # trace defaults for the display properties.
    warpByScalar1Display.ColorArrayName = ['CELLS', 'cell_level']
    warpByScalar1Display.LookupTable = celllevelLUT

    # hide data in view
    Hide(slice1, renderView1)
    UpdatePipeline()
    # show color bar/color legend
    warpByScalar1Display.SetScalarBarVisibility(renderView1, True)
    UpdatePipeline()
    # export view
    #ExportView('/home/milinda/Desktop/gr_sp1/slice_level_wbs.svg', view=renderView1)
    imgName=img_prefix+'_slice_level_wbs_'+str(imgIndex).zfill(6)+'.png'
    #ExportView(imgName, view=renderView1)
    #viewLayout1 = GetLayout()
    SaveScreenshot(imgName,magnification=1, quality=100,view=renderView1)
    #Delete(viewLayout1)
    #del viewLayout1
    print "image3 written"

    #slice1Display.SetScalarBarVisibility(renderView1, False)
    warpByScalar1Display.SetScalarBarVisibility(renderView1, False)
    # set scalar coloring
    ColorBy(warpByScalar1Display, ('POINTS', 'U_CHI'))

    # rescale color and/or opacity maps used to include current data range
    #warpByScalar1Display.RescaleTransferFunctionToDataRange(True)
    u_CHILUT = GetColorTransferFunction('U_CHI')
    u_CHILUT.LockDataRange = 0
    u_CHILUT.InterpretValuesAsCategories = 0
    u_CHILUT.ShowCategoricalColorsinDataRangeOnly = 0
    u_CHILUT.RescaleOnVisibilityChange = 0
    u_CHILUT.EnableOpacityMapping = 0
    u_CHILUT.RGBPoints = [6.286501118071329e-08, 0.231373, 0.298039, 0.752941, 0.4964692923160348, 0.865003, 0.865003, 0.865003, 0.9929385217670584, 0.705882, 0.0156863, 0.14902]
    u_CHILUT.UseLogScale = 0
    u_CHILUT.ColorSpace = 'Diverging'
    u_CHILUT.UseBelowRangeColor = 0
    u_CHILUT.BelowRangeColor = [0.0, 0.0, 0.0]
    u_CHILUT.UseAboveRangeColor = 0
    u_CHILUT.AboveRangeColor = [1.0, 1.0, 1.0]
    u_CHILUT.NanColor = [1.0, 1.0, 0.0]
    u_CHILUT.Discretize = 1
    u_CHILUT.NumberOfTableValues = 1024
    u_CHILUT.ScalarRangeInitialized = 1.0
    u_CHILUT.HSVWrap = 0
    u_CHILUT.VectorComponent = 0
    u_CHILUT.VectorMode = 'Magnitude'
    u_CHILUT.AllowDuplicateScalars = 1
    u_CHILUT.Annotations = []
    u_CHILUT.ActiveAnnotatedValues = []
    u_CHILUT.IndexedColors = []

    # toggle 3D widget visibility (only when running from the GUI)
    #Hide3DWidgets(proxy=slice1.SliceType)

    # get color legend/bar for u_CHILUT in view renderView1
    u_CHILUTColorBar = GetScalarBar(u_CHILUT, renderView1)
    u_CHILUTColorBar.AutoOrient = 1
    u_CHILUTColorBar.Orientation = 'Vertical'
    u_CHILUTColorBar.WindowLocation = 'LowerRightCorner'
    u_CHILUTColorBar.Position = [0.89, 0.02]
    u_CHILUTColorBar.Title = "U_CHI"
    u_CHILUTColorBar.ComponentTitle = ''
    u_CHILUTColorBar.TitleJustification = 'Centered'
    u_CHILUTColorBar.TitleColor = [1.0, 1.0, 1.0]
    u_CHILUTColorBar.TitleOpacity = 1.0
    u_CHILUTColorBar.TitleFontFamily = 'Arial'
    u_CHILUTColorBar.TitleBold = 0
    u_CHILUTColorBar.TitleItalic = 0
    u_CHILUTColorBar.TitleShadow = 0
    u_CHILUTColorBar.TitleFontSize = 16
    u_CHILUTColorBar.LabelColor = [1.0, 1.0, 1.0]
    u_CHILUTColorBar.LabelOpacity = 1.0
    u_CHILUTColorBar.LabelFontFamily = 'Arial'
    u_CHILUTColorBar.LabelBold = 0
    u_CHILUTColorBar.LabelItalic = 0
    u_CHILUTColorBar.LabelShadow = 0
    u_CHILUTColorBar.LabelFontSize = 16
    u_CHILUTColorBar.AutomaticLabelFormat = 0
    u_CHILUTColorBar.LabelFormat = '%-#2.2g'
    u_CHILUTColorBar.DrawTickMarks = 1
    u_CHILUTColorBar.DrawTickLabels = 1
    u_CHILUTColorBar.UseCustomLabels = 0
    u_CHILUTColorBar.CustomLabels = []
    u_CHILUTColorBar.AddRangeLabels = 1
    u_CHILUTColorBar.RangeLabelFormat = '%-#6.1e'
    u_CHILUTColorBar.DrawAnnotations = 1
    u_CHILUTColorBar.AddRangeAnnotations = 0
    u_CHILUTColorBar.AutomaticAnnotations = 0
    u_CHILUTColorBar.DrawNanAnnotation = 0
    u_CHILUTColorBar.NanAnnotation = 'NaN'
    u_CHILUTColorBar.TextPosition = 'Ticks right/top, annotations left/bottom'
    u_CHILUTColorBar.ScalarBarThickness = 16
    u_CHILUTColorBar.ScalarBarLength = 0.33

    # show color bar/color legend
    warpByScalar1Display.SetScalarBarVisibility(renderView1, True)
    # Rescale transfer function
    u_CHILUT.RescaleTransferFunction(0.0001, 1.0)
    u_CHIPWF.RescaleTransferFunction(0.0001, 1.0)
    UpdatePipeline()
    renderView1.Update()
    # export view
    #ExportView('/home/milinda/Desktop/gr_sp1/slice_wbs.svg', view=renderView1)
    imgName=img_prefix+'_slice_wbs_'+str(imgIndex).zfill(6)+'.png'
    #ExportView(imgName, view=renderView1)
    #viewLayout1 = GetLayout()
    SaveScreenshot(imgName,magnification=1, quality=100,view=renderView1)
    #Delete(viewLayout1)
    #del viewLayout1
    print "image4 written"
    # hide data in view
    #Hide(warpByScalar1, renderView1)

    # set active source
    #SetActiveSource(slice1)

    #Delete(uCHILUT)
    #del uCHILUT

    #Delete(uCHIPWF)
    #del uCHIPWF

    #Delete(celllevelLUT)
    #del celllevelLUT

    #Delete(celllevelPWF)
    #del celllevelPWF

    #Delete(bssn_gr_0pvtuDisplay)
    #del bssn_gr_0pvtuDisplay

    #Delete(slice1Display)
    #del slice1Display

    #Delete(warpByScalar1Display)
    #del warpByScalar1Display

    #Delete(annotateTimeFilter1Display)
    #del annotateTimeFilter1Display
    #SetActiveSource(annotateTimeFilter1)
    #Delete(annotateTimeFilter1)
    #del annotateTimeFilter1

    # destroy warpByScalar1
    SetActiveSource(warpByScalar1)
    Delete(warpByScalar1)
    del warpByScalar1

    # set active source
    #SetActiveSource(bssn_gr_0pvtu)

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

    #tracker.print_diff()
    #mem_end=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    #print 'Memory usage for current iteration: %s (kb)' % (mem_end-mem_begin)
    #print 'Memory usage use total (before gc): %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    #gc.collect()
    #print 'Memory usage use total (after gc): %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
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
    #step=0
    saveImg(file_prefix,img_prefix,step,imgCount)
    imgCount=imgCount+1;

#
