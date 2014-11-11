 # -*- coding: utf-8 -*-
"""
Create visualization with standard vtk actors, renders, windowsn, interactors

USAGE: 
=============
from display import *
loadDisplay = Display()  
loadDisplay.dicomTransform(image, image_pos_pat, image_ori_pat)
loadDisplay.addSegment(lesion3D)
loadDisplay.subImage(Images2Sub, timep)                  
loadDisplay.visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Methods:
=============
dicomTransform(image, image_pos_pat, image_ori_pat)
addSegment(lesion3D)
subImage(Images2Sub, timep)                  
visualize(images, image_pos_pat, image_ori_pat, sub, postS, interact)

Class Instance Attributes:
===============
'origin': (-167.0, -69.0, -145.0)
'spacing': (0.44920000433921814, 0.44920000433921814, 3.0)
'dims': (512, 512, 96), 

VTK Instance objects:
=============
'xImagePlaneWidget': (vtkImagePlaneWidget)
'yImagePlaneWidget': (vtkImagePlaneWidget)
'zImagePlaneWidget': (vtkImagePlaneWidget)
'picker': (vtkCellPicker)
'iren1': (vtkWin32RenderWindowInteractor)
'camera': (vtkOpenGLCamera)
'mapper_mesh': (vtkPainterPolyDataMapper)
'actor_mesh': (vtkOpenGLActor)
'renWin1': (vtkWin32OpenGLRenderWindow)
'renderer1': (vtkOpenGLRenderer)


Created on Tue Apr 01 10:18:34 2014
@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
from numpy import *
import re
import collections
import math

class Display(object):
    """
    USAGE:
    =============
    loadDisplay = Display()
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.picker = vtk.vtkCellPicker()
        self.picker.SetTolerance(0.005) 
        
        # Create 3 orthogonal view using the ImagePlaneWidget
        self.xImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.yImagePlaneWidget = vtk.vtkImagePlaneWidget()
        self.zImagePlaneWidget = vtk.vtkImagePlaneWidget()
        
        #  The 3 image plane widgets
        self.xImagePlaneWidget.DisplayTextOn();
        self.xImagePlaneWidget.SetPicker(self.picker);
        self.xImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.xImagePlaneWidget.SetKeyPressActivationValue('x');
        self.xImagePlaneWidget.GetPlaneProperty().SetColor(1, 0, 0);
        self.xImagePlaneWidget.SetResliceInterpolateToNearestNeighbour();
        
        self.yImagePlaneWidget.DisplayTextOn();
        self.yImagePlaneWidget.SetPicker(self.picker);
        self.yImagePlaneWidget.RestrictPlaneToVolumeOn();
        self.yImagePlaneWidget.SetKeyPressActivationValue('y');
        self.yImagePlaneWidget.GetPlaneProperty().SetColor(0, 1, 0);
        self.yImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        
        self.zImagePlaneWidget.DisplayTextOn();
        self.zImagePlaneWidget.SetPicker(self.picker);
        self.zImagePlaneWidget.SetKeyPressActivationValue('z');
        self.zImagePlaneWidget.GetPlaneProperty().SetColor(0, 0, 1);
        self.zImagePlaneWidget.SetLookupTable(self.xImagePlaneWidget.GetLookupTable());
        self.zImagePlaneWidget.SetRightButtonAutoModifier(1);
        
        # Create a renderer, render window, and render window interactor to
        # display the results.
        self.renderer1 = vtk.vtkRenderer()
        self.renWin1 = vtk.vtkRenderWindow()
        self.iren1 = vtk.vtkRenderWindowInteractor()
        
        self.renWin1.SetSize(1000, 800);
        self.renWin1.AddRenderer(self.renderer1)
        self.iren1.SetRenderWindow(self.renWin1)
        
        self.xImagePlaneWidget.SetInteractor( self.iren1 )
        self.yImagePlaneWidget.SetInteractor( self.iren1 )
        self.zImagePlaneWidget.SetInteractor( self.iren1 )
        
        # Set Up Camera view
        self.camera = self.renderer1.GetActiveCamera()
        self.renderer1.SetBackground(0.0, 0.0, 0.0)
        self.iren1.SetPicker(self.picker)
                
        self.T1origin = [0,0,0]
        self.T2origin = [0,0,0]
        self.T2extent = [0,0,0,0,0,0]
        self.T1extent = [0,0,0,0,0,0]
        self.T1spacing = [0,0,0]
        
    def __call__(self):       
        """ Turn Class into a callable object """
        Display()
        
        
    def addSegment(self, lesion3D, color, interact):        
        '''Add segmentation to current display'''
        # Set the planes based on seg bounds
        self.lesion_bounds = lesion3D.GetBounds()
        print "\n Mesh DICOM bounds: "
        print "xmin, xmax= [%d, %d]" % (self.lesion_bounds[0], self.lesion_bounds[1])
        print "yin, ymax= [%d, %d]" %  (self.lesion_bounds[2], self.lesion_bounds[3]) 
        print "zmin, zmax= [%d, %d]" % (self.lesion_bounds[4], self.lesion_bounds[5])
        
        ### GEt semgnetation information
        self.no_pts_segm = lesion3D.GetNumberOfPoints()
        print "no pts %d" % self.no_pts_segm
        
        # get VOI volume
        VOI_massProperty = vtk.vtkMassProperties()
        VOI_massProperty.SetInput(lesion3D)
        VOI_massProperty.Update()
               
        # VTK is unitless. The units you get out are the units you put in.
        # If your input polydata has points defined in terms of millimetres, then
        # the volume will be in cubic millimetres. 
        self.VOI_vol = VOI_massProperty.GetVolume() # mm3
        self.VOI_surface = VOI_massProperty.GetSurfaceArea() # mm2
    
        # just print the results
        print "\nVolume lesion = ", self.VOI_vol
        print "Surface lesion  = ", self.VOI_surface
        
        # Calculate the effective diameter of the surface D=2(sqrt3(3V/(4pi))) 
        diam_root = (3*self.VOI_vol)/(4*pi)
        self.VOI_efect_diameter = 2*pow(diam_root,1.0/3) 
        print "VOI_efect_diameter = ", self.VOI_efect_diameter
            
        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInput( lesion3D )
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()
        
        # centroid of lesion 
        self.lesion_centroid = [0,0,0]
        self.lesion_centroid = centerOfMassFilter.GetCenter()
        print "lesion_centroid = ", self.lesion_centroid
        self.lesioncentroidijk()
        
        # Add ICPinit_mesh.vtk to the render
        self.mapper_mesh = vtk.vtkPolyDataMapper()
        self.mapper_mesh.SetInput( lesion3D )
        self.mapper_mesh.ScalarVisibilityOff()
        
        self.actor_mesh = vtk.vtkActor()
        self.actor_mesh.SetMapper(self.mapper_mesh)
        self.actor_mesh.GetProperty().SetColor(color)    #R,G,B
        self.actor_mesh.GetProperty().SetOpacity(0.3)
        self.actor_mesh.GetProperty().SetPointSize(5.0)
        self.actor_mesh.GetProperty().SetRepresentationToWireframe()
        
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetSliceIndex( 0 )
        
        self.renderer1.AddActor(self.actor_mesh)
        
        # Initizalize
        self.renderer1.Modified()
        self.renWin1.Render()
        self.renderer1.Render()
        
        if(interact==True):
            self.iren1.Start()
                
        return 
        
        
    def subImage(self, Images2Sub, timep):
        '''subtract volumes based on indicated postS'''
        sub_preMat = vtk.vtkImageMathematics()
        sub_preMat.SetOperationToSubtract()
        sub_preMat.SetInput1(Images2Sub[timep])
        sub_preMat.SetInput2(Images2Sub[0])
        sub_preMat.Update()
                    
        sub_pre = vtk.vtkImageData()
        sub_pre =  sub_preMat.GetOutput()
        # define image based on subtraction of postS -preS
        subtractedImage = sub_pre
        
        return subtractedImage
      
      
    def visualize(self, ImagedataVTK, sub, postS, interact):
        '''Display and render volumes, reference frames, actors and widgets'''
        if(sub):
            #subtract volumes based on indicated postS            
            # define image based on subtraction of postS -preS
            self.transformed_image = self.subImage(images, postS)
        else:
            self.transformed_image = ImagedataVTK[postS]# images[postS]            
             
        # Proceed to build reference frame for display objects based on DICOM coords   
        #[self.transformed_image, transform_cube] = self.dicomTransform(image, image_pos_pat, image_ori_pat)
        # get info from image before visualization
        self.transformed_image.UpdateInformation()
        self.dims = self.transformed_image.GetDimensions()
        print "Image Dimensions"
        print self.dims
        self.T1spacing = self.transformed_image.GetSpacing()
        print "Image Spacing"
        print self.T1spacing
        self.T1origin = self.transformed_image.GetOrigin()
        print "Image Origin"
        print self.T1origin
        self.T1extent = list(self.transformed_image.GetWholeExtent())
        print "Image Extent"
        print self.T1extent
            
        # Set up ortogonal planes
        self.xImagePlaneWidget.SetInput( self.transformed_image )
        self.xImagePlaneWidget.SetPlaneOrientationToXAxes()
        self.xImagePlaneWidget.SetSliceIndex(0)
        self.yImagePlaneWidget.SetInput( self.transformed_image )
        self.yImagePlaneWidget.SetPlaneOrientationToYAxes()
        self.yImagePlaneWidget.SetSliceIndex(0)
        self.zImagePlaneWidget.SetInput( self.transformed_image )
        self.zImagePlaneWidget.SetPlaneOrientationToZAxes()
        self.zImagePlaneWidget.SetSliceIndex(0)
            
        self.xImagePlaneWidget.On()
        self.yImagePlaneWidget.On()
        self.zImagePlaneWidget.On()
        
        ############
        # bounds and initialize camera
        bounds = self.transformed_image.GetBounds()
        self.renderer1.ResetCamera(bounds)    
        self.renderer1.ResetCameraClippingRange()
        self.camera.SetViewUp(0.0,-1.0,0.0)
        self.camera.Azimuth(315)
        
        # Initizalize
        self.renWin1.Modified()
        self.renWin1.Render()
        self.renderer1.Render()
        
        if(interact==True):
            interactor = self.renWin1.GetInteractor()
            interactor.Start()
                            
        return
        
    def lesioncentroidijk(self):
        im_pt=[0,0,0]
        ijk=[0,0,0]
        pco = [0,0,0]
        pixId_sliceloc = self.transformed_image.FindPoint(self.lesion_centroid)
        self.transformed_image.GetPoint(pixId_sliceloc, im_pt) 
        io = self.transformed_image.ComputeStructuredCoordinates( im_pt, ijk, pco)
        if io:
            self.lesion_centroid_ijk = ijk
            print "\n Lesion centroid"
            print self.lesion_centroid_ijk
            
        return
         
        
    def extract_segment_dims(self, lesion3D):
        '''Extract mayor dimensions of automatic lesion segmentation for validation'''
        # define image based on subtraction of postS -preS          
        axis_lenghts = array( [0,0,0] ).astype(float)
        
        l_bounds = lesion3D.GetBounds()
        axis_lenghts[0] = (l_bounds[1]-l_bounds[0])**2     
        axis_lenghts[1] = (l_bounds[3]-l_bounds[2])**2
        axis_lenghts[2] = (l_bounds[5]-l_bounds[4])**2 
        
        return axis_lenghts
        