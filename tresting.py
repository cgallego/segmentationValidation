# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 13:18:32 2014

@author: windows
"""
from vtk.util.numpy_support import vtk_to_numpy
import pickle
import vtk

def display(vtkimage, interact):
    picker = vtk.vtkCellPicker()
    picker.SetTolerance(0.005) 
    
    # Create 3 orthogonal view using the ImagePlaneWidget
    xImagePlaneWidget = vtk.vtkImagePlaneWidget()
    yImagePlaneWidget = vtk.vtkImagePlaneWidget()
    zImagePlaneWidget = vtk.vtkImagePlaneWidget()
    
    #  The 3 image plane widgets
    xImagePlaneWidget.DisplayTextOn();
    xImagePlaneWidget.SetPicker(picker);
    xImagePlaneWidget.RestrictPlaneToVolumeOn();
    xImagePlaneWidget.SetKeyPressActivationValue('x');
    xImagePlaneWidget.GetPlaneProperty().SetColor(1, 0, 0);
    xImagePlaneWidget.SetResliceInterpolateToNearestNeighbour();
    
    yImagePlaneWidget.DisplayTextOn();
    yImagePlaneWidget.SetPicker(picker);
    yImagePlaneWidget.RestrictPlaneToVolumeOn();
    yImagePlaneWidget.SetKeyPressActivationValue('y');
    yImagePlaneWidget.GetPlaneProperty().SetColor(0, 1, 0);
    yImagePlaneWidget.SetLookupTable(xImagePlaneWidget.GetLookupTable());
    
    zImagePlaneWidget.DisplayTextOn();
    zImagePlaneWidget.SetPicker(picker);
    zImagePlaneWidget.SetKeyPressActivationValue('z');
    zImagePlaneWidget.GetPlaneProperty().SetColor(0, 0, 1);
    zImagePlaneWidget.SetLookupTable(xImagePlaneWidget.GetLookupTable());
    zImagePlaneWidget.SetRightButtonAutoModifier(1);
    
    # Create a renderer, render window, and render window interactor to
    # display the results.
    renderer1 = vtk.vtkRenderer()
    renWin1 = vtk.vtkRenderWindow()
    iren1 = vtk.vtkRenderWindowInteractor()
    
    renWin1.SetSize(1000, 800);
    renWin1.AddRenderer(renderer1)
    iren1.SetRenderWindow(renWin1)
    
    xImagePlaneWidget.SetInteractor( iren1 )
    yImagePlaneWidget.SetInteractor( iren1 )
    zImagePlaneWidget.SetInteractor( iren1 )
    
    # Set Up Camera view
    camera = renderer1.GetActiveCamera()
    renderer1.SetBackground(0.0, 0.0, 0.0)
    iren1.SetPicker(picker)
    
    dims = vtkimage.GetDimensions()
    print "Image Dimensions"
    print dims
    T1spacing = vtkimage.GetSpacing()
    print "Image Spacing"
    print T1spacing
    T1origin = vtkimage.GetOrigin()
    print "Image Origin"
    print T1origin
    T1extent = list(vtkimage.GetWholeExtent())
    print "Image Extent"
    print T1extent
        
    # Set up ortogonal planes
    xImagePlaneWidget.SetInput( vtkimage )
    xImagePlaneWidget.SetPlaneOrientationToXAxes()
    xImagePlaneWidget.SetSliceIndex(0)
    yImagePlaneWidget.SetInput( vtkimage )
    yImagePlaneWidget.SetPlaneOrientationToYAxes()
    yImagePlaneWidget.SetSliceIndex(0)
    zImagePlaneWidget.SetInput( vtkimage )
    zImagePlaneWidget.SetPlaneOrientationToZAxes()
    zImagePlaneWidget.SetSliceIndex(0)
        
    xImagePlaneWidget.On()
    yImagePlaneWidget.On()
    zImagePlaneWidget.On()
    
    ############
    # bounds and initialize camera
    bounds = vtkimage.GetBounds()
    renderer1.ResetCamera(bounds)    
    renderer1.ResetCameraClippingRange()
    camera.SetViewUp(0.0,-1.0,0.0)
    camera.Azimuth(315)
       
    # Initizalize    
    if(interact==True):
        # Initialize the interactor.
        iren1.Initialize()
        renWin1.Render()
        iren1.Start()
    
    return renderer1, renWin1, iren1
    
    
def convertSegmentation2vtkImage(npMask, templateImage):
    """ Takes a numpy.ndarray and converts it to a vtkimageData. require npImagesandMask to pass on image info """
    # Create vtk object
    #dims=shape(npMask)
    dims = templateImage.GetDimensions()
    size_array = dims[0]*dims[1]*dims[2]
    flatim = npMask.flatten()
    
    # create vtk image
    vtk_image = vtk.vtkImageData()
    vtk_image.DeepCopy(templateImage)
    vtk_image.SetNumberOfScalarComponents(1)
    vtk_image.SetScalarTypeToDouble()
    vtk_image.AllocateScalars()
    
    # create vtk double array
    image_array = vtk.vtkDoubleArray() 
    image_array.SetNumberOfComponents(1) 
    
    # not too efficient convertion of np.array to vtk. Far from ideal
    for k in range(size_array):
        image_array.InsertTuple1(k,float(flatim[k]))
        
    vtk_image.GetPointData().SetScalars(image_array) 
    vtk_image.Update()
      
    return vtk_image

def vtkimage2poly(vtkIm):
    thresh = vtk.vtkImageThreshold()
    thresh.SetInput(vtkIm)
    thresh.ThresholdByUpper(0.5)
    thresh.SetInValue(255.0)
    thresh.SetOutValue(0.0)
    thresh.Update()
      
    contouriso = vtk.vtkMarchingCubes()
    contouriso.SetInput(thresh.GetOutput())
    contouriso.SetValue(0,255)
    contouriso.ComputeScalarsOn()
    contouriso.SetNumberOfContours(1)
    contouriso.Update()
    
    # Recalculate num_voxels and vol_lesion on VOI
    segmpoly = contouriso.GetOutput()
    nvoxels = segmpoly.GetNumberOfCells()
    npoints = segmpoly.GetNumberOfPoints()
    print "Number of points: %d" % npoints 
    print "Number of cells: %d" % nvoxels 
        
    return segmpoly
    
 
def addSegment(segmpoly, color, interact, renderer1, renWin1, iren1):        
    '''Add segmentation to current display'''
    # Set the planes based on seg bounds
    segmpoly_bounds = segmpoly.GetBounds()
    print "\n Mesh DICOM bounds: "
    print "xmin, xmax= [%d, %d]" % (segmpoly_bounds[0], segmpoly_bounds[1])
    print "yin, ymax= [%d, %d]" %  (segmpoly_bounds[2], segmpoly_bounds[3]) 
    print "zmin, zmax= [%d, %d]" % (segmpoly_bounds[4], segmpoly_bounds[5])
    
    ### GEt semgnetation information
    no_pts_segm = segmpoly.GetNumberOfPoints()
    print "no pts %d" % no_pts_segm
    
    # get VOI volume
    VOI_massProperty = vtk.vtkMassProperties()
    VOI_massProperty.SetInput(segmpoly)
    VOI_massProperty.Update()
           
    # VTK is unitless. The units you get out are the units you put in.
    # If your input polydata has points defined in terms of millimetres, then
    # the volume will be in cubic millimetres. 
    VOI_vol = VOI_massProperty.GetVolume() # mm3
    VOI_surface = VOI_massProperty.GetSurfaceArea() # mm2

    # just print the results
    print "\nVolume lesion = ", VOI_vol
    print "Surface lesion  = ", VOI_surface
    
    # Calculate the effective diameter of the surface D=2(sqrt3(3V/(4pi))) 
    diam_root = (3*VOI_vol)/(4*pi)
    VOI_efect_diameter = 2*pow(diam_root,1.0/3) 
    print "VOI_efect_diameter = ", VOI_efect_diameter
        
    centerOfMassFilter = vtk.vtkCenterOfMass()
    centerOfMassFilter.SetInput( segmpoly )
    centerOfMassFilter.SetUseScalarsAsWeights(False)
    centerOfMassFilter.Update()
    
    # centroid of lesion 
    lesion_centroid = [0,0,0]
    lesion_centroid = centerOfMassFilter.GetCenter()
    print "lesion_centroid = ", lesion_centroid
            
    # Add ICPinit_mesh.vtk to the render
    mapper_mesh = vtk.vtkPolyDataMapper()
    mapper_mesh.SetInput( segmpoly )
    mapper_mesh.ScalarVisibilityOff()
    
    actor_mesh = vtk.vtkActor()
    actor_mesh.SetMapper(mapper_mesh)
    actor_mesh.GetProperty().SetColor(color)    #R,G,B
    actor_mesh.GetProperty().SetOpacity(0.3)
    actor_mesh.GetProperty().SetPointSize(5.0)
    actor_mesh.GetProperty().SetRepresentationToWireframe()
    
    renderer1.AddActor(actor_mesh)
    
    if(interact==True):
        iren1.Initialize()
        renWin1.Render()
        iren1.Start()
            
    return 
    
    
#templateImage =vtk.vtkImageData()
#templateImage.SetDimensions([512,512,48])
#templateImage.SetSpacing([1,1,3])
#dims = templateImage.GetDimensions()
##### for test
#myarr = zeros((dims))
#myarr[240:280, 240:280, 15:25]=1
#myarrt = myarr.transpose(2,1,0)
#myimvtkimage = convertSegmentation2vtkImage(myarrt, templateImage)


#### forreal
templateImage=load.DICOMImages[0]
npMasktvtkimage = convertSegmentation2vtkImage(npMask, templateImage)

# marching cubes on image
segmpoly = vtkimage2poly(npMasktvtkimage)
[renderer1, renWin1, iren1] = display(npMasktvtkimage, interact=False)
addSegment(segmpoly, [0,1,0], True, renderer1, renWin1, iren1)