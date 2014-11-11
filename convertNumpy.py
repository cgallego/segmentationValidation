# -*- coding: utf-8 -*-
"""
Convert numpy data into vtk objects

USAGE:
=============
from convertNumpy import *
convertImages = convertNumpy()
convertImages.convertArray2vtkImage( nparray, npImagesandMask )
convertImages.convertDCEArray2vtkImage(npImagesandMask, postS)


Class Methods:
=============
convertArray2vtkImage(npImagesandMask, postS)
createMeshfromMask(npmask, npImagesandMask)
convertDCEArray2vtkImage(npImagesandMask, postS)

Class Instance Attributes:
===============
'vtk_image' vtkImageData
'meshlesion3D'  vtkPolyData
'npDICOMImages' dictionary: when running testVTK2pOutputfile()

dict_keys([
'im0', 	float64 array of pre-contrast volume
'im1', 	float64 array of post-contrast1 volume
'im2',	float64 array of post-contrast2 volume 
'im3', 	float64 array of post-contrast3 volume
'im4', 	float64 array of post-contrast4 volume
'mask', float64 array 3D volume of lesion segmentation, zeros in background, ones in foreground 

'ti0' 	dicom tag [0x0008,0x0032] for Volume Series im0
'ti1', 	dicom tag [0x0008,0x0032] for Volume Series im1
'ti2', 	dicom tag [0x0008,0x0032] for Volume Series im2
'ti3', 	dicom tag [0x0008,0x0032] for Volume Series im3
'ti4', 	dicom tag [0x0008,0x0032] for Volume Series im4
'image_ori_pat', 	dicom tag [0x0020,0x0037]
'image_pos_pat', 	dicom tag [0x0020,0x0032] from most-far-left slice
'spacing', tuple 3	[inplane, slice_thickness]
'dims',    tuple 3	[Volume extent]
'nvol',    int	Number of total dynamic volumes (e.g 5)
])


Created on Fri Apr 04 11:36:56 2014
@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
from numpy import *
from display import *
from inputs_init import *
from vtk.util.numpy_support import vtk_to_numpy
import pickle

#!/usr/bin/env python
class convertNumpy(object):
    """
    USAGE:
    =============
    convertImages = convertNumpy()
    convertImages.convertArray2vtkImage(npImagesandMask, postS)
    """
    def __init__(self): 
        """ Initialize image container """
        self.meshlesion3D = []
                
    def __call__(self):       
        """ Turn Class into a callable object """
        convertNumpy()
        
    def convertSegmentation2vtkImage(self, npMask, templateImage):
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
        
        
    def createMeshfromMask(self, vtkimMask):
        """ Takes a binary mask image and makes it a vtkPolyData VOI_mesh """                
        # Create a binary Image with 0-255
        image_VOIlesion = vtk.vtkImageThreshold()
        image_VOIlesion.ThresholdByUpper(0.5)
        image_VOIlesion.SetInValue(255)
        image_VOIlesion.SetOutValue(0)
        image_VOIlesion.SetInput(vtkimMask)
        image_VOIlesion.Update()
        
        # Convert VOIlesion into polygonal struct
        VOIlesion_poly = vtk.vtkMarchingCubes() 
        VOIlesion_poly.SetValue(0,125)
        VOIlesion_poly.SetInput(image_VOIlesion.GetOutput())
        VOIlesion_poly.ComputeNormalsOff()
        VOIlesion_poly.Update()
        
        # prepare output 
        self.meshlesion3D = VOIlesion_poly.GetOutput()
        
        # Recalculate num_voxels and vol_lesion on VOI
        nvoxels = self.meshlesion3D.GetNumberOfCells()
        npoints = self.meshlesion3D.GetNumberOfPoints()
        print "Number of points: %d" % npoints 
        print "Number of cells: %d" % nvoxels 
                
        return self.meshlesion3D
        
        
    def convertDCEArray2numpy(self, DICOMimage):
        
        dims = DICOMimage.GetDimensions()
        spacing = DICOMimage.GetSpacing()
        im_scalars = DICOMimage.GetPointData().GetScalars()
        np_imdata = vtk_to_numpy(im_scalars) 
        np_imdata = np_imdata.reshape(dims[2], dims[1], dims[0]) 
        #np_imdata = array(np_imdata.transpose(2,1,0)).astype(float) 
        # append
        numpy_image = np_imdata
    
        return numpy_image, dims, spacing
        
        
    def convertDCEArray2vtkImage(self, numpy_image, dims, spacing, side):        
        #####################
        # Create vtk object
        print numpy_image.shape
        size_array = dims[0]*dims[1]*dims[2]
        #numpy_image = numpy_image[::-1,::-1,:]
        flatim = numpy_image.flatten()
        flatim = flatim[::-1]
        
        # create vtk image
        vtk_image = vtk.vtkImageData()
        vtk_image.SetSpacing(spacing)
        vtk_image.SetDimensions(dims[0],dims[1],dims[2])
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
        center = vtk_image.GetCenter()
        print center
        
        # deal with bottom-left order
        flipx = vtk.vtkImageFlip()
        flipx.SetInput(vtk_image)  
        flipx.SetFilteredAxis(0)
        flipx.FlipAboutOriginOff()
        flipx.ReleaseDataFlagOn()
        flipx.Update()            
        
        flipz = vtk.vtkImageFlip()
        flipz.SetInput(flipx.GetOutput())  
        flipz.SetFilteredAxis(2)
        flipz.FlipAboutOriginOff()
        flipz.ReleaseDataFlagOn()
        flipz.Update() 
        
        vtk_image = flipz.GetOutput()
        print vtk_image.GetOrigin()
        
        return vtk_image
        
        
    def testconvertArray2vtkImage(self, dims):        
        #####################
        # Create vtk object
        narray = zeros(dims)
        narray[5:15,5:15,4:7]=1 # create a white rectangle in the middle
        size_array = narray.size
        flatim = narray.transpose(2,1,0)
        flatim = flatim.reshape(dims[2], dims[1], dims[0]) 
        flatim = flatim.flatten()
        
        # create vtk image
        vtk_image = vtk.vtkImageData()
        vtk_image.SetDimensions(20,20,10)
        vtk_image.SetOrigin(0,0,0)
        vtk_image.SetSpacing(1,1,2)
        vtk_image.SetScalarTypeToUnsignedShort()
        vtk_image.AllocateScalars()
        # create vtk double array
        image_array = vtk.vtkDoubleArray() 
        image_array.SetNumberOfComponents(1) 
        
        # not too efficient convertion of np.array to vtk. Far from ideal
        for k in range(size_array):
            image_array.InsertNextTupleValue([flatim[k]])
            
        vtk_image.GetPointData().SetScalars(image_array) 
        vtk_image.Update()
             
        
        return vtk_image
    
    def testVTK2pOutputfile(self):
        """ Takes a Study/ExamID/SeriesID/ dicom images and convertes them into a data structure to test integration pipeline"""
        # Open filename list
        StudyID = '18'    
        DicomExamNumber = '7714' # corresponds to old way of ret
        Lesions_id = '1721'
        SeriesID = 'S44' # corresponds to dynamic sequence;
            
        ###### Loading 
        print "Start by loading volumes..."
        load = Inputs_init()
        [series_path, phases_series, lesionID_path] = load.readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
        print "Path to series location: %s" % series_path 
        print "List of pre and post contrast volume names: %s" % phases_series
        print "Path to lesion segmentation: %s" % lesionID_path
        
        print "\n Load Segmentation..."
        lesion3D = load.loadSegmentation(lesionID_path)
        print "Data Structure: %s" % lesion3D.GetClassName()
        print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
        print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
        
        print "\n Visualize volumes..."
        loadDisplay = Display()
        lesion3D_mesh = loadDisplay.addSegment(lesion3D)
        loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=3, interact=False)

        #######################################################
        ###### Testing integration format change of input data 
        #######################################################    
        # Convert load.DICOMImages data to list of arrays [x,y,z] and lesion3D segmentation to mask [x,y,z]
        self.npDICOMImages = {}
        for i in range(len(load.DICOMImages)):
            # convert 'DICOMImages': list[(vtkImageData) to npDICOMImages': list[(ndarray)
            dims = load.DICOMImages[i].GetDimensions()
            spacing = load.DICOMImages[i].GetSpacing()
            im_scalars = load.DICOMImages[i].GetPointData().GetScalars()
            np_imdata = vtk_to_numpy(im_scalars) 
            np_imdata = np_imdata.reshape(dims[2], dims[1], dims[0]) 
            np_imdata = array(np_imdata.transpose(2,1,0)).astype(float) 
            # append
            self.npDICOMImages['im'+str(i)] = np_imdata
            
            # process time points needed for dynamic features
            abspath_PhaseID = series_path+os.sep+str(phases_series[i]) 
            # Get total number of files
            [len_listSeries_files, FileNms_slices_sorted_stack] = processDicoms.ReadDicomfiles(abspath_PhaseID)
            mostleft_slice = FileNms_slices_sorted_stack.slices[0]
            
            # Get dicom header, retrieve
            dicomInfo_series = dicom.read_file(abspath_PhaseID+os.sep+str(mostleft_slice))
            # (0008,0032) AT S Acquisition Time       # hh.mm.ss.frac
            ti = str(dicomInfo_series[0x0008,0x0032].value)
            self.npDICOMImages['ti'+str(i)]=ti
            
            
        # create other information from dicom data
        self.npDICOMImages['dims'] = load.DICOMImages[0].GetDimensions()
        self.npDICOMImages['spacing'] = load.DICOMImages[0].GetSpacing()
        self.npDICOMImages['nvol'] = len(load.DICOMImages)
        self.npDICOMImages['image_pos_pat'] = load.image_pos_pat # position of far most left (indicates origin)
        self.npDICOMImages['image_ori_pat'] = load.image_ori_pat
        
        ################################################################ NEEDED TO TEST CHANGING FORMAT OF DATA
        # Create mask for VOI
        [transformed_image, t] = Display().dicomTransform(load.DICOMImages[0], load.image_pos_pat, load.image_ori_pat)
        self.vtkmask = load.createVTKMaskfromMesh(lesion3D, transformed_image) # SHOULD RETURN A VTKIMAGEDATA REPRESENTING MASK
        
        # save image as metafile image
        vtkimage_w = vtk.vtkMetaImageWriter()
        vtkimage_w.SetInput(transformed_image)
        vtkimage_w.SetFileName( 'vtkimage.mhd' )
        vtkimage_w.Write()
        
        # ## save mask as metafile image
        vtkmask_w = vtk.vtkMetaImageWriter()
        vtkmask_w.SetInput(self.vtkmask )
        vtkmask_w.SetFileName( 'vtkmask.mhd' )
        vtkmask_w.Write()
        
        # write to image        
        maskscalars = self.vtkmask.GetPointData().GetScalars()
        npmask = vtk_to_numpy(maskscalars)     
        npmask = npmask.reshape(self.npDICOMImages['dims'][2], self.npDICOMImages['dims'][1], self.npDICOMImages['dims'][0]) 
        npmask = array(npmask.transpose(2,1,0)).astype(float) 
                
        self.npDICOMImages['mask'] = npmask  # SHOULD RETURN A NUMPY ARRAY REPRESENTING MASK
        
        # Save a dictionary into a pickle file. to retrieve later
        # Not saving the arrays corectly
        pickle.dump( self.npDICOMImages, open( "npDICOMImages.p", "wb" ), -1 )
        
        ###################################################### FINISH TESTING
    
        return 
    
