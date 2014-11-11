# -*- coding: utf-8 -*-
"""
 USAGE:
=============
from inputs_init import *
input = Inputs_init()
input.readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
input.loadSegmentation(lesionID_path)

Class Methods:
=============
readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
loadSegmentation(lesionID_path)

Class Instance Attributes:
===============
'slice_thickn': 3.0, 
'image_ori_pat': ['-0', '1', '0', '-0', '-0', '-1'], 
'spacing': (0.44920000433921814, 0.44920000433921814, 3.0), 
'dims': (512, 512, 96), 
'image_pos_pat': ['145.059', '-167.095', '69.3364'], 
'DICOMImages': list[(vtkImageData)
'lesion3D_mesh': (vtkPolyData)
    
Created on Wed Apr 02 13:42:50 2014

@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""

import os, os.path
import sys
import string
import dicom
import vtk

import processDicoms
import shlex, subprocess

#!/usr/bin/env python
class Inputs_init:
    """
    USAGE:
    =============
    input = Inputs_init()
    input.readVolumes(StudyID, DicomExamNumber, SeriesID, Lesions_id)
    
    """
    def __init__(self):
        self.DICOMImages = []
    
    
    def __call__(self):       
        """ Turn Class into a callable object """
        Inputs_init()
        
        # Gets only_files in directory of folder mydir, excluding subdirectories folders
    def get_only_filesindirectory(self, mydir):
        return [name for name in os.listdir(mydir) 
            if os.path.isfile(os.path.join(mydir, name))]
       
    
    def processVolumes(self, data_loc, StudyID, DicomExamNumber):
        # ask for which series to load 
        [abspath_ExamID, eID, SeriesIDall, studyFolder, dicomInfo] = processDicoms.get_series(StudyID, data_loc)
        
        print "\n----------------------------------------------------------"
        choseSerie = raw_input('Enter n DCE pre-contrast Series to load (0-n), or x to pass:')
        if (choseSerie != 'x'):
            ## append collection of cases    
            SeriesID = SeriesIDall[int(choseSerie)]
            print "\nSeriesID: %s" %  SeriesID
        else:
            SeriesID=""

        return SeriesID
        
        
    def readVolumes(self, data_loc, StudyID, DicomExamNumber, SeriesID, chooseSide, side):
        """
        ARGUMENTS:
        =============
        StudyID (str)           Study CAD patient number      
        DicomExamNumber (str)   Imaging Study Number    
        SeriesID  (str)         Dynamic Series number (e.g S600, S3)     
        
        OUTPUTS:
        =============
        series_path (str)       Path to series location      
        phases_series (list)      list of pre and post contrast volume names      
        lesionID_path  (str)      path to lesion segmentation      
        DICOMImages list(vtkImageData)      List of vtkImageData objects corresponding to dynamic series
        
        """
        # Get Root folder ( the directory of the script being run)
        series_path = data_loc+'/'+str(StudyID)+'/'+str(DicomExamNumber)
        print series_path
        # test
        ###############################################################
        os.chdir(series_path)   
        # Get series to load
        series_path = data_loc+os.sep+str(StudyID)+os.sep+str(DicomExamNumber)            
        phases_series=[]
        testSID = str(SeriesID)
        if 'S' in str(testSID):
            #print testSID[1:]
            chosen_phase = int(testSID[1:])
        else:
            chosen_phase = int(testSID)
        
        if(testSID[0] == 'S'):
            phases_series.append('S'+str(chosen_phase))
                            
            for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3, chosen_phase+4]:
                phases_series.append( 'S'+str(chSer) )    
        else:
            phases_series.append(str(chosen_phase))
                            
            for chSer in [chosen_phase+1, chosen_phase+2, chosen_phase+3, chosen_phase+4]:
                phases_series.append( str(chSer) )

        if os.path.exists('SeriesPhases'):
            print '''SeriesPhases'''

            # Get series to load
            series_path = data_loc+os.sep+str(StudyID)+os.sep+str(DicomExamNumber)+os.sep+str(SeriesID)
                                                        
            # process all Volumes when in stacks of Dyn Volumes
            if os.path.exists(series_path+os.sep+'pre-Contrast'):
                phases_series = []
                phases_series.append('pre-Contrast')
                        
                #"Arranging series scans"
                for i in range(1,5):
                    phases_series.append('post_Contrast-'+str(i))
        
        # Get total number of files and some DICOM tags needed fro DICOM coords
        if(chooseSide):
            # create left and right slip
            for i in range(0,len(phases_series)):
                abspath_PhaseID = series_path+os.sep+phases_series[i] 
                print abspath_PhaseID
                [len_listSeries_files, FileNms_slices_sorted_stack] = processDicoms.ReadDicomfiles(abspath_PhaseID)
                newlen_listSeries_files = int(round(len_listSeries_files/2))
                                 
                os.chdir(abspath_PhaseID)
                ## if AccessionN folder doesn't exist create it        
                if not os.path.exists("left"):
                    os.makedirs("left")
                if not os.path.exists("right"):
                    os.makedirs("right")
                    
                ## process left
                os.chdir("left")
                for k in range(newlen_listSeries_files):
                    cmd = 'mv ../'+str(FileNms_slices_sorted_stack.iloc[k]['slices'])+' '+str(FileNms_slices_sorted_stack.iloc[k]['slices'])
                    print "cmd -> " + cmd
                    p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
                    p1.wait()
                    
                os.chdir(abspath_PhaseID)
                os.chdir("right")
                for k in range(newlen_listSeries_files,len_listSeries_files):
                    cmd = 'mv ../'+str(FileNms_slices_sorted_stack.iloc[k]['slices'])+' '+str(FileNms_slices_sorted_stack.iloc[k]['slices'])
                    print "cmd -> " + cmd
                    p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE)
                    p1.wait()
                
            os.chdir(series_path)    
            pre_abspath_PhaseID = series_path+'/'+phases_series[0]+'/'+side
            [len_listSeries_files, FileNms_slices_sorted_stack] = processDicoms.ReadDicomfiles(pre_abspath_PhaseID)
            newlen_listSeries_files = int(round(len_listSeries_files/2))
            # read only side
            if(side == "left"):
                mostextreme_slice = FileNms_slices_sorted_stack.iloc[0]['slices']
                print FileNms_slices_sorted_stack.iloc[0]['location']
                # Get dicom header, retrieve: image_pos_pat and image_ori_pat
                dicomInfo_series=dicom.read_file(pre_abspath_PhaseID+os.sep+str(mostextreme_slice)) 
                
            if(side == "right"):
                mostextreme_slice = FileNms_slices_sorted_stack.iloc[newlen_listSeries_files]['slices']
                print FileNms_slices_sorted_stack.iloc[newlen_listSeries_files]['location']
                # Get dicom header, retrieve: image_pos_pat and image_ori_pat
                dicomInfo_series=dicom.read_file(pre_abspath_PhaseID+os.sep+str(mostextreme_slice)) 
            
            # get series info
            self.slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
            self.image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
            self.image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)
            
            for i in range(0,len(phases_series)):
                abspath_PhaseID = series_path+'/'+phases_series[i]+'/'+side
                print abspath_PhaseID
                                         
                # Done reading volume filenames now load them to DICOMimageReader                                         
                os.chdir(abspath_PhaseID)               
                dicomReader  = vtk.vtkDICOMImageReader()
                dicomReader.SetDirectoryName( abspath_PhaseID )
                dicomReader.Update()
                
                im = vtk.vtkImageData()
                im =  dicomReader.GetOutput()
                self.dims = im.GetDimensions()
                self.spacing = im.GetSpacing()
                
                print "VTK Dimensions im.GetDimensions(): %d %d %d" % self.dims
                print "VTK Spacing im.GetSpacing(): %f %f %f\n" % self.spacing
                
                # Append to objects image            
                self.DICOMImages.append( im )
                os.chdir(series_path)
                
        else:                    
            pre_abspath_PhaseID = series_path+'/'+phases_series[0]
            [len_listSeries_files, FileNms_slices_sorted_stack] = processDicoms.ReadDicomfiles(pre_abspath_PhaseID)
            mostleft_slice = FileNms_slices_sorted_stack.slices[0]
            
            # Get dicom header, retrieve: image_pos_pat and image_ori_pat
            dicomInfo_series = dicom.read_file(pre_abspath_PhaseID+os.sep+str(mostleft_slice)) 
            self.slice_thickn = float(dicomInfo_series[0x0018,0x0050].value)  
            
            # Image Position (0020,0032): specifies the x, y, and z coordinates of the upper left hand corner of the image. This tag specifies the coordinates 
            # of the the first voxel transmitted.
            self.image_pos_pat = list(dicomInfo_series[0x0020,0x0032].value)
            # Image Orientation (0020,0037): specifies the direction cosines 
            # of the first row and the first column with respect to the patient. 
            # The direction of the axes are defined by the patients orientation 
            # to ensure LPS system ( x-axis increasing to the left hand side of the patient, 
            # y-axis increasing to the posterior side of the patient and z-axis increasing toward
            # the head of the patient )
            self.image_ori_pat = list(dicomInfo_series[0x0020,0x0037].value)
            
            for i in range(0,len(phases_series)):
                abspath_PhaseID = series_path+os.sep+phases_series[i] 
                print abspath_PhaseID
                                         
                # Done reading volume filenames now load them to DICOMimageReader                                         
                os.chdir(abspath_PhaseID)               
                dicomReader  = vtk.vtkDICOMImageReader()
                dicomReader.SetDirectoryName( abspath_PhaseID )
                dicomReader.Update()
                
                im = vtk.vtkImageData()
                im =  dicomReader.GetOutput()
                self.dims = im.GetDimensions()
                self.spacing = im.GetSpacing()
                
                print "VTK Dimensions im.GetDimensions(): %d %d %d" % self.dims
                print "VTK Spacing im.GetSpacing(): %f %f %f\n" % self.spacing
                
                # Append to objects image            
                self.DICOMImages.append( im )
                os.chdir(data_loc)   
                        
        return(series_path, phases_series)        
        
    
    def loadSegmentation(self, lesionID_path):
        """
        ARGUMENTS:
        =============
        lesionID_path: (str)        path to lesion segmentation

        OUTPUT:
        =============
        lesion3D (vtkPolyData)      3D lesion segmentation as a vtkPolyData object       
        """ 
        # need to locate and read Lesion Seg
        VOIlesion = 'VOIlesion_selected.vtk'
            
        lesion3D_reader = vtk.vtkPolyDataReader()
        lesion3D_reader.SetFileName( lesionID_path+os.sep+VOIlesion )
        lesion3D_reader.Update()
        
        # Extract the polydata from the PolyDataReader        
        self.lesion3D_mesh = lesion3D_reader.GetOutput()
      
        return self.lesion3D_mesh
        
    def createVTKMaskfromMesh(self, VOI_mesh, im):
        """ Takes an image and a VOI_mesh and returns a boolean image with only 1s inside the VOI_mesh """                
        
        # Create an Image of Fext
        white_image = vtk.vtkImageData()
        white_image.DeepCopy(im) 
        white_image.SetScalarTypeToUnsignedChar()
        white_image.AllocateScalars()
        count = white_image.GetNumberOfPoints()
        for i in range(count):
            white_image.GetPointData().GetScalars().SetTuple1(i, 1.0)
  
        # polygonal data --> image stencil:
        pol2stenc = vtk.vtkPolyDataToImageStencil()
        pol2stenc.SetInput(VOI_mesh)
        pol2stenc.SetOutputOrigin(im.GetOrigin())
        pol2stenc.SetOutputSpacing(im.GetSpacing())
        pol2stenc.SetOutputWholeExtent(im.GetExtent())
        pol2stenc.Update()
         
        # cut the corresponding white image and set the background:
        imgstenc = vtk.vtkImageStencil()
        imgstenc.SetInput(white_image)
        imgstenc.SetStencil(pol2stenc.GetOutput())
        imgstenc.ReverseStencilOff()
        imgstenc.SetBackgroundValue(0.0)
        imgstenc.Update()
        
        return imgstenc.GetOutput()
        
    def createMeshfromVTKMask(self, imMask):
        """ Takes a binary mask image and makes it a vtkPolyData VOI_mesh """                
        
        # Create a binary Image with 0-255
        image_VOIlesion = vtk.vtkImageThreshold()
        image_VOIlesion.ThresholdByUpper(0.1)
        image_VOIlesion.SetInValue(255)
        image_VOIlesion.SetOutValue(0)
        image_VOIlesion.SetInput(imMask)
        image_VOIlesion.Update()
        
        # Convert VOIlesion into polygonal struct
        VOIlesion_poly = vtk.vtkMarchingCubes() 
        VOIlesion_poly.SetValue(0,125)
        VOIlesion_poly.SetInput(image_VOIlesion.GetOutput())
        VOIlesion_poly.ComputeNormalsOff()
        VOIlesion_poly.Update()
        
        # Recalculate num_voxels and vol_lesion on VOI
        nvoxels = VOIlesion_poly.GetOutput().GetNumberOfCells()
        npoints = VOIlesion_poly.GetOutput().GetNumberOfPoints()
        print "Number of points: %d" % npoints 
        print "Number of cells: %d" % nvoxels 
        
        # prepare output 
        meshlesion3D = VOIlesion_poly.GetOutput()
        
        return meshlesion3D
        
       
    
    
    
    
    
    
    
    
    
    
    
    
        
        
    