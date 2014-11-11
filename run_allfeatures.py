# -*- coding: utf-8 -*-
"""
Created on Tue Apr 01 10:47:40 2014

This script will load volumes, Load a lesion Segmentation (VOI), Visualize volumes
and then extract Dynamic, Morphology and Texture features from the VOI.

@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""
from inputs_init import *
from dictionaries import data_loc
import pickle

from convertNumpy import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
from add_records import *

###################################################### 
if __name__ == '__main__':    
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__)) 
    print path_rootFolder
    lesion_id=40
    
    # Open filename list
    file_ids = open(sys.argv[1],"r")
    for fileline in file_ids:
        # Get the line: Study#, DicomExam#
        fileline = fileline.split()
        PatientID = "CADpat" 
        cond = fileline[0]
        StudyID = fileline[1]
        AccessionN = fileline[2]
        SeriesID = fileline[3]
        massornm = fileline[5]
        segid = fileline[4]
        side = fileline[6]
        lesion_id +=1
        print "\n Lesion %d " % lesion_id
        
        # Load the dictionary back from the pickle file.
        print "\n Loading Data..."
        pth_pickle = 'Z:/Hongbo/segmentedlesions/'+cond+'-'+massornm+'-VOIlesion_id'+StudyID+'_'+AccessionN+'_'+SeriesID+'_'+segid+'.mask'
        pkl_file = open(pth_pickle, 'rb')
        npMask = pickle.load(pkl_file)
        pkl_file.close()
        
        load = Inputs_init()
        #SeriesID = load.processVolumes(data_loc, StudyID, AccessionN)
        ###### Start by Loading 
        print "Start by loading volumes..."
        [series_path, phases_series] = load.readVolumes(data_loc, StudyID, AccessionN, SeriesID, chooseSide=True, side=side)
        print "Path to series location: %s" % series_path 
        print "List of pre and post contrast volume names: %s" % phases_series
        
        # at this point proceed with all image data needed stored in npMask and mask data in meshlesion3D    
        # Convert back to vtk objects for actual use
         # create convertor instance
        convertImages = convertNumpy()
        templateImage=load.DICOMImages[0]
        ImagedataNumpy=[] 
        ImagedataVTK=[]
        for i in range( len(load.DICOMImages) ):
            [numpy_image, dims, spacing] = convertImages.convertDCEArray2numpy(load.DICOMImages[i])
            ImagedataNumpy.append( numpy_image )
            # Convert back to vtk objects for actual use
            vtk_image = convertImages.convertDCEArray2vtkImage(numpy_image, dims, spacing, side)
            ImagedataVTK.append(vtk_image)
             
        # Convert segmentation
        npMasktvtk_image = convertImages.convertSegmentation2vtkImage(npMask, templateImage)

        ## creat VOI 3D mesh from binary mask
        print "\n Load Segmentation..."
        meshlesion3D = convertImages.createMeshfromMask(npMasktvtk_image)
        
        # transform to dicom to get same coords for mask     
        print "\n Reload and visualize"
        loadDisplay = Display()
        loadDisplay.visualize(ImagedataVTK, sub=False, postS=1, interact=False)
        loadDisplay.addSegment(meshlesion3D, (0,1,0), interact=True)
        
        #############################
        ###### 4) Extract Dynamic features
        #############################
        loadDynamic = Dynamic()
        print "\n Extract Dynamic contour features..."
        dyn_contour = loadDynamic.extractfeatures_contour(ImagedataVTK, side, series_path, phases_series, meshlesion3D)
        print "\n=========================================="
        print dyn_contour

        print "\n Extract Dynamic inside features..."
        dyn_inside = loadDynamic.extractfeatures_inside(ImagedataVTK, side, series_path, phases_series, meshlesion3D)
        print dyn_inside
        print "\n=========================================="
         
        #############################
        ###### 5) Extract Morphology features
        #############################
        loadMorphology = Morphology()
        print "\n Extract Morphology features..."
        morphofeatures = loadMorphology.extractfeatures(ImagedataVTK, series_path, phases_series, meshlesion3D)
        print "\n=========================================="
        print morphofeatures
        print "\n=========================================="
        
        #############################        
        ###### 6) Extract Texture features
        #############################
        loadTexture = Texture()
        print "\n Extract Texture features..."
        texturefeatures = loadTexture.extractfeatures(ImagedataVTK, series_path, phases_series, meshlesion3D, loadMorphology.VOI_efect_diameter, loadMorphology.lesion_centroid )
        print "\n=========================================="
        print texturefeatures
        print "\n=========================================="
        
        #############################
        ###### Send record to DB
        ## append collection of cases
        #############################  
        records = AddRecords()
        if(cond=="malignant"):  label=massornm+'M'
        if(cond=="benign"):  label=massornm+'B'
        
        label
        print "\n Adding record case to DB..."
        records.lesion_2DB(segid, StudyID, 'NA', AccessionN, datetime.date(9999, 12, 31), 'NA', 
                           'NA', 1, 0, side, 'NA', 
                            datetime.date(9999, 12, 31), 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', label,  cond)
                            
        if "mass" == massornm:
            records.mass_2DB(lesion_id, cond, SeriesID, 'NA', 'NA', 'NA' )

        if "nonmass" == massornm: 
            records.nonmass_2DB(lesion_id, cond, SeriesID, 'NA', 'NA', 'NA')
    
        # send features
        # Dynamic
        records.dyn_records_2DB(lesion_id, dyn_inside['A.inside'], dyn_inside['alpha.inside'], dyn_inside['beta.inside'], dyn_inside['iAUC1.inside'], dyn_inside['Slope_ini.inside'], dyn_inside['Tpeak.inside'], dyn_inside['Kpeak.inside'], dyn_inside['SER.inside'], dyn_inside['maxCr.inside'], dyn_inside['peakCr.inside'], dyn_inside['UptakeRate.inside'], dyn_inside['washoutRate.inside'], dyn_inside['maxVr.inside'], dyn_inside['peakVr.inside'], dyn_inside['Vr_increasingRate.inside'], dyn_inside['Vr_decreasingRate.inside'], dyn_inside['Vr_post_1.inside'],
                               dyn_contour['A.contour'], dyn_contour['alpha.contour'], dyn_contour['beta.contour'], dyn_contour['iAUC1.contour'], dyn_contour['Slope_ini.contour'], dyn_contour['Tpeak.contour'], dyn_contour['Kpeak.contour'], dyn_contour['SER.contour'], dyn_contour['maxCr.contour'], dyn_contour['peakCr.contour'], dyn_contour['UptakeRate.contour'], dyn_contour['washoutRate.contour'], dyn_contour['maxVr.contour'], dyn_contour['peakVr.contour'], dyn_contour['Vr_increasingRate.contour'], dyn_contour['Vr_decreasingRate.contour'], dyn_contour['Vr_post_1.contour'] )
        
        # Morphology
        records.morpho_records_2DB(lesion_id, morphofeatures['min_F_r_i'], morphofeatures['max_F_r_i'], morphofeatures['mean_F_r_i'], morphofeatures['var_F_r_i'], morphofeatures['skew_F_r_i'], morphofeatures['kurt_F_r_i'], morphofeatures['iMax_Variance_uptake'], 
                                                  morphofeatures['iiMin_change_Variance_uptake'], morphofeatures['iiiMax_Margin_Gradient'], morphofeatures['k_Max_Margin_Grad'], morphofeatures['ivVariance'], morphofeatures['circularity'], morphofeatures['irregularity'], morphofeatures['edge_sharp_mean'],
                                                  morphofeatures['edge_sharp_std'], morphofeatures['max_RGH_mean'], morphofeatures['max_RGH_mean_k'], morphofeatures['max_RGH_var'], morphofeatures['max_RGH_var_k'] )
        # Texture
        records.texture_records_2DB(lesion_id, texturefeatures['texture_contrast_zero'], texturefeatures['texture_contrast_quarterRad'], texturefeatures['texture_contrast_halfRad'], texturefeatures['texture_contrast_threeQuaRad'], 
                                                  texturefeatures['texture_homogeneity_zero'], texturefeatures['texture_homogeneity_quarterRad'], texturefeatures['texture_homogeneity_halfRad'], texturefeatures['texture_homogeneity_threeQuaRad'], 
                                                  texturefeatures['texture_dissimilarity_zero'], texturefeatures['texture_dissimilarity_quarterRad'], texturefeatures['texture_dissimilarity_halfRad'], texturefeatures['texture_dissimilarity_threeQuaRad'], 
                                                  texturefeatures['texture_correlation_zero'], texturefeatures['texture_correlation_quarterRad'], texturefeatures['texture_correlation_halfRad'], texturefeatures['texture_correlation_threeQuaRad'], 
                                                  texturefeatures['texture_ASM_zero'], texturefeatures['texture_ASM_quarterRad'], texturefeatures['texture_ASM_halfRad'], texturefeatures['texture_ASM_threeQuaRad'], 
                                                  texturefeatures['texture_energy_zero'], texturefeatures['texture_energy_quarterRad'], texturefeatures['texture_energy_halfRad'], texturefeatures['texture_energy_threeQuaRad'] )
        
        # SEgmentation details
        records.segment_records_2DB(lesion_id, loadDisplay.lesion_bounds[0], loadDisplay.lesion_bounds[1], loadDisplay.lesion_bounds[2], loadDisplay.lesion_bounds[3], loadDisplay.lesion_bounds[4], loadDisplay.lesion_bounds[5],
                                    loadDisplay.no_pts_segm, loadDisplay.VOI_vol, loadDisplay.VOI_surface, loadDisplay.VOI_efect_diameter, str(list(loadDisplay.lesion_centroid)), str(list(loadDisplay.lesion_centroid_ijk)))
                                                    
    file_ids.close()
    
    
    
    
    
    
        
        
    

