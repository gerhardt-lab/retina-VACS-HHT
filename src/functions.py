"""
functions for the loading, analysing and saving data (might be split into separate files later)
import functions
"""

# imports (TODO: remove unnecessary)
from __future__ import division
from skimage import io # change import style
from pylab import np, plt
from skimage.io import imread, imshow, imsave
#import io
from skimage.filters import rank, gaussian, threshold_otsu
from skimage.morphology import disk, remove_small_objects, skeletonize, label
import scipy.misc
import scipy
import seaborn as sns
import pandas as pd
import os
from fnmatch import fnmatch
import csv
import time
import scipy.misc
from shutil import copyfile
from scipy.stats import binned_statistic
#from joblib import Parallel, delayed
import multiprocessing
import shutil
import yaml
import random
import matplotlib.colors as colors

sns.set_context('poster')
sns.set_style('white')


def get_parameters():

    with open('conf/base/parameters.yml') as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        parameter_dict = yaml.load(file, Loader=yaml.FullLoader)
       
    with open('conf/local/parameters.yml') as file:
        parameter_list_local = yaml.load(file, Loader=yaml.FullLoader)

        for key in parameter_list_local:
            parameter_dict[key] = parameter_list_local[key]
        
        print("Parameters: ")
        print(parameter_dict)
        
    return parameter_dict


    
def get_key_file(parameters):
    
    key_file = pd.read_excel(parameters['key_file'], skiprows=[0])
        
    # adapt filename / file paths according to the os system
        
    for i, row in key_file.iterrows():
        path  = str(row["filename"]).replace("\\ ", " ")
        path  = path.replace("\#", "#")
        key_file.at[i,"filename"] = path
        
    # key_file.to_csv("data/processed_key_file.csv")
            
    return key_file

def get_tiff_meta_data(file_path):
    
    import PIL 
    
    meta_data = dict()
    
    float_attributes = ["PhysicalSizeX","PhysicalSizeY"]
    int_attributes = ["SizeX","SizeY"]
    
    img = PIL.Image.open(file_path)
    
    #print(img.tag.tagdata)
    for key in img.tag.tagdata:
   
        
        tag_data = str(img.tag.tagdata[key])
        
        #print(entry)
        if tag_data.find("xml")>0:
            #print(entry)
            #print(entry.decode("utf-8"))
            
        #    print(img.tag.tagdata[key])
            data = tag_data.split(" ")
            #print(data)
            
            #print(entries)
            for d in data:
                attribute = d.split("=")[0]
                if attribute in float_attributes:
                    meta_data[attribute] = float(d.split("=")[1].replace('"',''))
                if attribute in int_attributes:
                    meta_data[attribute] = int(d.split("=")[1].replace('"',''))
                                            
                    
    return meta_data

def TP_from_filename(fn):
    if 'P5_P6' in fn:
        return 'P5_P6'
    elif 'P5_P7' in fn:
        return 'P5_P7'
    elif 'P5_P8' in fn:
        return 'P5_P8'
    elif 'P5_P9' in fn:
        return 'P5_P9'
    elif 'P5_P11' in fn:
        return 'P5_P11'
    elif 'P5_P12' in fn:
        return 'P5_P12'
    elif 'P8_P9' in fn:
        return 'P8_P9'
    elif 'P3_P7' in fn:
        return 'P3_P7'
    elif 'P3_P15' in fn:
        return 'P3_P15'
    elif 'P5_P15' in fn:
        return 'P5_P15'
    else:
        print('TP_from_filename: Could not extract known time point from ' + str(fn) + '.')
        return None
    
def assign_color(tp):
    if tp == 'P5_P6':
        return 'orange'
    if tp == 'P5_P7':
        return 'blue'
    if tp == 'P5_P8':
        return 'gray'
    if tp == 'P5_P9':
        return 'green'
    if tp == 'P5_P11':
        return'red'
    if tp == 'P5_P12':
        return 'yellow'
    if tp == 'P8_P9':
        return 'purple'
    if tp == 'P3_P7':
        return 'purple'
    else:
        print('assign_color: Could not extract known time point from ' + str(tp) + '.')
        return None

def expID_from_results_path(fn):
    try:
        folderName = fn.split("python_results")[1]
        expID = int(folderName.split('_')[1])
        return expID
    except (IndexError, ValueError):
        print('Experiment ID could not be extracted from path:')
        print(fn)

def get_image_type_from_path(fn):
    try:
        folderName = fn.split("python_results")[1]
        if 'GFP X ERG/' in folderName:
            return 'GFP X ERG'
        elif 'GFP/' in folderName:
            return 'GFP'
        elif '/ERG/' in folderName:
            return 'ERG'
        else:
            print('Microscopy file type (GFP, ERG or GFP X ERG) could not be extracted from path:')
            print(fn)
            return None
    except (IndexError, ValueError):
        print('Microscopy file type (GFP, ERG or ERG X GFP) could not be extracted from path:')
        print(fn)
        return None

def get_injectionDay(tp):
    try:
        return int(tp.split('_')[0].split('P')[1])
    except:
        print('get_injectionDay: Could not extract time point from ' + str(tp) + '.')
        return None

def get_cullingDay(tp):
    try:
        return int(tp.split('_')[1].split('P')[1])
    except:
        print('get_cullingDay: Could not extract known time point from ' + str(tp) + '.')
        return None

def get_migrationTime(tp):
    try:
        tp1 = tp.split("_")[0]
        tp2 = tp.split("_")[1]
        return int(tp2.split('P')[1])-int(tp1.split('P')[1])
    except:
        print('get_migrationTime: Could not extract time point from ' + str(tp) + '.')
        return None

def get_radial_indices(radii, minimum, maximum):
    if minimum == 0 and maximum != 0:
        return np.where(radii<maximum)[0]
    elif minimum == 0 and maximum == 0:
        return np.where(radii==maximum)[0]
    else:
        indices_minimum = np.where(radii<=minimum)[0]
        indices_maximum = np.where(radii<maximum)[0]
        #return np.array(list(filter(lambda x: x not in indices_minimum, indices_maximum)))
        return list(set(indices_maximum.flat)-set(indices_minimum.flat))

#def save_all_positions(parameters, key_file, experimentID, drawn_filename, gfp_filename, gfpXerg_filename, erg_filename):
def save_all_positions(parameters, key_file, experimentID, drawn_filename, microscopy_filename, attributes, dir_append):
    #img_map = np.array(io.imread(folder_name + '/drawn.ome.tif')).astype('bool')
    #bin_img = io.imread(folder_name + '/GFP_threshold.tif').astype('bool')
    
    folder_name = parameters["out_dir"] + 'processed/python_results/'
    print(microscopy_filename)
    # get time point to include in output folder later. This helps when reading in the results in load_res
    #tp = TP_from_filename(drawn_filename)
    experiment_df = key_file[key_file["ExperimentID"] == experimentID]
    tp = 'P' + str(experiment_df["THO Injection Point"].iloc[0]) + '_P' + str(experiment_df["Collection Point"].iloc[0])
    
    dir_name = folder_name + 'ExperimentID_' + str(experimentID) + '_' + str(tp) + '/' + dir_append

    # if directory exists already, don't compute anything
    if os.path.exists(dir_name):
        print('Results for experiment ID %s (%s) were computed already. Moving to next file.' %(experimentID, dir_append))
    else:
        # print info of experiment that is used for computing results
        experiment_df  = key_file[key_file["ExperimentID"] == experimentID]
        print(experiment_df[experiment_df["Drawn"] == 1].iloc[0])
    
        img_map = np.array(io.imread(drawn_filename)).astype('bool')
        print(img_map.shape)
        # problem that accoured with Yi's retinas when they were RGB
        if img_map.ndim != 3:
            print('WARNING! Something was wrong with the draw file dimensions...')
            img_map = np.moveaxis(img_map[:,:,:,0], 0, -1)
        # masks that included the sprounting front showed a different dimension
        if img_map.shape[0] < img_map.shape[2]:
            print('WARNING! Something was wrong with the draw file dimensions.')
            print('Old dimensions: ' + str(img_map.shape))
            img_map = np.moveaxis(img_map[:,:,:], 0, -1)
            print('New dimensions: ' + str(img_map.shape))
        if img_map.shape[2] > 4:
            if np.count_nonzero(img_map[:,:,4] == True) > (img_map.shape[0]*img_map.shape[1]/2):
                img_map[:,:,4] = np.invert(img_map[:,:,4])
            img_map[:,:,4] = np.invert(scipy.ndimage.morphology.binary_fill_holes(img_map[:,:,4]))

        # print mask just to make sure everything's alright
        print('Do these masks look okay?')
        fig, axes = plt.subplots(1, img_map.shape[2], figsize=(10, 8.2), sharex=True, sharey=True)
        
        for i in range(img_map.shape[2]):
            axes[i].imshow(img_map[:,:,i])
            axes[i].axis('off')  
        plt.close()
        
        # make sure that the microscopy image is not inverted (i.e. more than half the pixels are True)
        bin_img = np.array(io.imread(microscopy_filename)).astype('bool')
        if bin_img.ndim == 3:
            print('Warning! The microscopy image was loaded with 3 channels instead of just 1. :(')
            bin_img = bin_img[:,:,0]
        if np.count_nonzero(bin_img == True) > (bin_img.shape[0]*bin_img.shape[1]/2):
            bin_img = np.invert(bin_img)
        print(bin_img.shape)
        bin_img_flat = bin_img.flatten()
             
        
        # im img_map, the underlying structures are saved in the following way:
        pos_on = 0  # optical nerve
        pos_artery = 1
        pos_vein = 2
        pos_mask = 3
        pos_sf = 4 # sprouting front

    
        # define distance map of distances to veins, arteries and optical nerve, and account for pixel size
        # if um/pixel ratio cannot be extracted from any TIFF, take it from the key file instead
        #if not ("PhysicalSizeX" in attributes):
        key_single_exp = key_file[key_file["ExperimentID"] == experimentID]
        attributes["PhysicalSizeX"] = key_single_exp["Pixel size in um"].iloc[0]
            #print("Warning: scaling pixels to microns could not be extracted from TIFF and was taken from key file!!!")
        print("Attributes")
        print(attributes)
    
        scale_factor = attributes["PhysicalSizeX"]
        
        artery_distance = scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_artery]))*scale_factor
        artery_distance_flat = artery_distance.flatten()
        vein_distance = scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_vein]))*scale_factor
        vein_distance_flat = vein_distance.flatten()
        optical_distance = scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[:, :, pos_on]))*scale_factor
        optical_distance_flat = optical_distance.flatten()

        # define arterial-veinal distance
        av_distance = vein_distance / (artery_distance + vein_distance)
        av_distance[np.isnan(av_distance)] = 0.5
        av_distance_flat = av_distance.flatten()
    
        # define vectors for binning (edges and midpoints)
        v_AD = np.linspace(0, 3500, 51)
        v_AD_mps = (v_AD[1:] + v_AD[:-1]) / 2.
        v_VD = np.linspace(0, 3500, 51)
        v_VD_mps = (v_VD[1:] + v_VD[:-1]) / 2.
        v_R = np.linspace(0, 3500, 51)
        v_R_mps = (v_R[1:] + v_R[:-1]) / 2.
        v_AV = np.linspace(0, 1, 31)
        v_AV_mps = (v_AV[1:] + v_AV[:-1]) / 2.
    
        # calculate 1D histograms and means for GFP distribution
        distances_AD = artery_distance_flat[bin_img_flat]
        distances_VD = vein_distance_flat[bin_img_flat]
        distances_R = optical_distance_flat[bin_img_flat]
        distances_AV = av_distance_flat[bin_img_flat]
    
        hist_AD = np.histogram(distances_AD,
                                   bins=v_AD)
        mean_AD = np.mean(distances_AD)
    
        hist_VD = np.histogram(distances_VD,
                                   bins=v_VD)
        mean_VD = np.mean(distances_VD)
    
        hist_R = np.histogram(distances_R,
                                  bins=v_R)
        mean_R = np.mean(distances_R)
    
        hist_AV = np.histogram(distances_AV,
                                   bins=v_AV)
        mean_AV = np.mean(distances_AV)
    
        # calculate 2D histograms for GFP distribution
        hist_AD_R = np.histogram2d(optical_distance_flat[bin_img_flat],
                                       artery_distance_flat[bin_img_flat],
                                       bins=[v_R, v_AD])
        hist_VD_R = np.histogram2d(optical_distance_flat[bin_img_flat],
                                       vein_distance_flat[bin_img_flat],
                                       bins=[v_R, v_VD])
        hist_AV_R = np.histogram2d(optical_distance_flat[bin_img_flat],
                                       av_distance_flat[bin_img_flat],
                                       bins=[v_R, v_AV])
    
        # save all the results in .npy files in a separate directory
        os.makedirs(dir_name)
        
        np.save(dir_name + '/hist_AD', hist_AD)
        np.save(dir_name + '/hist_VD', hist_VD)
        np.save(dir_name + '/hist_AV', hist_AV)
        np.save(dir_name + '/hist_R', hist_R)
        np.save(dir_name + '/hist_AD_R', hist_AD_R)
        np.save(dir_name + '/hist_VD_R', hist_VD_R)
        np.save(dir_name + '/hist_AV_R', hist_AV_R)
        np.save(dir_name + '/hist_mps_AD', v_AD_mps)
        np.save(dir_name + '/hist_mps_VD', v_VD_mps)
        np.save(dir_name + '/hist_mps_AV', v_AV_mps)
        np.save(dir_name + '/hist_mps_R', v_R_mps)
        np.save(dir_name + '/mean_AV', mean_AV)
        np.save(dir_name + '/mean_VD', mean_VD)
        np.save(dir_name + '/mean_AD', mean_AD)
        np.save(dir_name + '/mean_R', mean_R)
    
        np.save(dir_name + '/distances_AD', distances_AD)
        np.save(dir_name + '/distances_VD', distances_VD)
        np.save(dir_name + '/distances_AV', distances_AV)
        np.save(dir_name + '/distances_R', distances_R)
        
        if img_map.shape[2] > 4:
            sproutingFront_distance = scipy.ndimage.morphology.distance_transform_edt(np.invert(img_map[ :,:, pos_sf]))
            sproutingFront_distance_flat = sproutingFront_distance.flatten()
            v_SF = np.linspace(0, 3500, 51)
            v_SF_mps = (v_SF[1:] + v_SF[:-1]) / 2.
            distances_SF = sproutingFront_distance_flat[bin_img_flat]
            np.save(dir_name + '/distances_SF', distances_SF)
        
        # save results in return struct as well ->
        ret_struct = {'retina_name': dir_name,
                      'v_AD_mps': v_AD_mps, 'v_VD_mps': v_VD_mps, 'v_AV_mps': v_AV_mps, 'v_R_mps': v_R_mps,
                      'hist_AD': hist_AD, 'hist_VD': hist_VD, 'hist_AV': hist_AV,
                      'hist_R': hist_R,
                      'hist_AD_R': hist_AD_R, 'hist_VD_R': hist_VD_R, 'hist_AV_R': hist_AV_R,
                      'mean_AD': mean_AD, 'mean_VD': mean_VD, 'mean_AV': mean_AV,
                      'mean_R': mean_R}
    
        print("Job done: " + dir_name)


        print("Plotting masks and distance maps to check...")
        folder = parameters['out_dir'] + "distance_maps/"
        if not os.path.exists(folder):
            os.makedirs(folder)
        # plot distance transforms for checking
        outline_mask = scipy.ndimage.morphology.binary_fill_holes(img_map[:, :, pos_mask])
        fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)
        cmap1 = colors.ListedColormap(['none', 'white'])
        axes[0].imshow(np.ma.masked_where(outline_mask == 0, optical_distance), cmap=plt.get_cmap("Greys"))
        # axes[0].imshow(optical_distance, cmap=plt.get_cmap("Greys"))
        axes[0].axis('off')
        CS = axes[0].contour(optical_distance, [1000, 1500, 2000], linewidths=2)
        axes[0].clabel(CS, [1000, 1500, 2000], inline=1, fmt='%1.1f', fontsize=12)
        axes[0].imshow(img_map[:, :, 0], cmap=cmap1)

        axes[1].imshow(np.ma.masked_where(outline_mask == 0, artery_distance), cmap=plt.get_cmap("Reds_r"))
        # axes[1].imshow(artery_distance, cmap=plt.get_cmap("Reds_r"))
        axes[1].axis('off')
        CS = axes[1].contour(optical_distance, [1000, 1500, 2000], linewidths=2)
        axes[1].clabel(CS, [1000, 1500, 2000], inline=1, fmt='%1.1f', fontsize=12)
        axes[1].imshow(img_map[:, :, 1], cmap=cmap1)

        axes[2].imshow(np.ma.masked_where(outline_mask == 0, vein_distance), cmap=plt.get_cmap("Blues_r"))
        # axes[2].imshow(vein_distance, cmap=plt.get_cmap("Blues_r"))
        axes[2].axis('off')
        CS = axes[2].contour(optical_distance, [1000, 1500, 2000], linewidths=2)
        axes[2].clabel(CS, [1000, 1500, 2000], inline=1, fmt='%1.1f', fontsize=12)
        axes[2].imshow(img_map[:, :, 2], cmap=cmap1)

        palette = plt.cm.coolwarm
        palette.set_bad('w', 1.0)
        axes[3].imshow(np.ma.masked_where(outline_mask == 0, av_distance), interpolation='bilinear', cmap=palette,
                       norm=colors.Normalize(vmin=0.0, vmax=1.0))
        # axes[3].imshow(av_distance, interpolation='bilinear', cmap=palette,
        #              norm=colors.Normalize(vmin=0.0, vmax=1.0))
        axes[3].axis('off')
        CS = axes[3].contour(optical_distance, [1000, 1500, 2000], linewidths=2)
        axes[3].clabel(CS, [1000, 1500, 2000], inline=1, fmt='%1.1f', fontsize=12)
        axes[3].imshow(img_map[:, :, 1], cmap=cmap1)

        plt.savefig(folder + 'distance_maps_ExpID-' + str(experimentID) + '.png', format="png", bbox_inches="tight",
                    dpi=150)
        plt.clf()
        plt.close('all')
        del fig, axes

        # plot masks for checking
        # find center of mass of optic nerve
        center_on = scipy.ndimage.measurements.center_of_mass(
            scipy.ndimage.morphology.binary_fill_holes(img_map[:, :, pos_on]))
        folder =  parameters['out_dir'] + "check_masks_dir/"
        if not os.path.exists(folder):
            os.makedirs(folder)
        fig, axes = plt.subplots(1, img_map.shape[2], figsize=(10, 8.2), sharex=True, sharey=True)

        for i in range(img_map.shape[2]):
            circle = plt.Circle((center_on[0], center_on[1]), 1500 / scale_factor, fill=False, linestyle='--')
            axes[i].imshow(img_map[:, :, i])
            axes[i].scatter(center_on[0], center_on[1], s=30, marker='+')
            axes[i].add_patch(circle)
            axes[i].axis('off')
        plt.savefig(folder + 'mask_ExpID-' + str(experimentID) + '.png', format="png", bbox_inches="tight", dpi=150)
        plt.clf()
        plt.close('all')

        # plot GFP image
        folder =  parameters['out_dir'] + "check_GFP_image_dir/"
        if not os.path.exists(folder):
            os.makedirs(folder)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.imshow(bin_img)
        plt.savefig(folder + 'microscopy-image_ExpID-' + str(experimentID) + '.png', format="png",
                    bbox_inches="tight", dpi=150)
        plt.clf()
        plt.close('all')

        return ret_struct



def compute_res(parameters, key_file):
    """ 
    go through all directories in input_dir_list, find all subfolders,
    compute the area, erg and GFP distributions and save the results there in .npy files
    
    input:
        - parameters
    return:

    """
    
    for experimentID in key_file["ExperimentID"].unique():
        
        experiment_df  = key_file[key_file["ExperimentID"] == experimentID]
        
        row_drawn = experiment_df[experiment_df["Drawn"] == 1].iloc[0]
        row_gfp = experiment_df[experiment_df["GFP"] == 1].iloc[0]

        drawn_filename = parameters["data_dir"] + row_drawn['filename']
        gfp_filename = parameters["data_dir"] + row_gfp['filename']

        # there are not always GFPxERG and ERG files, so we need to check and set the file names to False if there are none
        #if 1 in experiment_df["GFP X ERG"].values:
         #   row_gfpXerg = experiment_df[experiment_df["GFP X ERG"] == 1].iloc[0]
          #  gfpXerg_filename = parameters["data_dir"] + row_gfpXerg['filename']
        #else:
        gfpXerg_filename = False

        #if 1 in experiment_df["ERG"].values:
         #   row_erg = experiment_df[experiment_df["ERG"] == 1].iloc[0]
          #  erg_filename = parameters["data_dir"] + row_erg['filename']
        #else:
        erg_filename = False

        attributes = {}
        if gfp_filename != False:
            attributes = get_tiff_meta_data(gfp_filename)
        if erg_filename != False:
            if not ("PhysicalSizeX" in attributes):
                attributes = get_tiff_meta_data(erg_filename)
        if gfpXerg_filename != False:
            if not ("PhysicalSizeX" in attributes):
                attributes = get_tiff_meta_data(gfpXerg_filename)

        #print(row_drawn)
        if gfp_filename != False:
            save_all_positions(parameters, key_file, experimentID, drawn_filename, gfp_filename, attributes, 'GFP')
        if erg_filename != False:
            save_all_positions(parameters, key_file, experimentID, drawn_filename, erg_filename, attributes, 'ERG')
        if gfpXerg_filename != False:
            save_all_positions(parameters, key_file, experimentID, drawn_filename, gfpXerg_filename, attributes, 'GFP X ERG')


def load_res(parameters, dir_append):

    """ load all the results into on dataframe after they have been computed using compute_res
    
    input:
        - parameters
        - dir_append: string, either "GFP", "GFP X ERG" or "ERG" depending on
          which results are to be loaded (loading everything results in a dataframe that is too big)
    return:
        - df : pandas data frame
    """

    #inputDir = "data/processed/python_results/"
    inputDir = parameters["out_dir"] + "processed/python_results/" 

    folderList = []
    for path, subdirs, files in os.walk(inputDir):
        if 'ExperimentID' in path and dir_append in path:
            folderList.append(path)
    folder_list = np.unique(folderList)
    
    """
    there is a "ValueError: Object arrays cannot be loaded when allow_pickle=False" with the normal np.load function.
    I found this workaround. the np.load function is restored after the for loop
    """
    # save np.load
    #np_load_old = np.load
    np.load.__defaults__=(None, True, True, 'ASCII')
    
    # modify the default parameters of np.load
    #np.load = lambda *a, **k: np_load_old(*a, allow_pickle=True)#, **k)

    res_list = []
    for res_folder in folder_list:
        res = {'retina_name': res_folder + '/',
               'v_AD_mps': np.load(res_folder + '/hist_mps_AD.npy'),
               'v_VD_mps': np.load(res_folder + '/hist_mps_VD.npy'),
               'v_AV_mps': np.load(res_folder + '/hist_mps_AV.npy'),
               'v_R_mps': np.load(res_folder + '/hist_mps_R.npy'),
               'hist_AD_GFP': np.load(res_folder + '/hist_AD.npy'),
               'hist_VD_GFP': np.load(res_folder + '/hist_VD.npy'),
               'hist_AV_GFP': np.load(res_folder + '/hist_AV.npy'),
               'hist_R_GFP': np.load(res_folder + '/hist_R.npy'),
               'hist_AD_R_GFP': np.load(res_folder + '/hist_AD_R.npy'),
               'hist_VD_R_GFP': np.load(res_folder + '/hist_VD_R.npy'),
               'hist_AV_R_GFP': np.load(res_folder + '/hist_AV_R.npy'),
               'mean_AD_GFP': np.load(res_folder + '/mean_AD.npy'),
               'mean_VD_GFP': np.load(res_folder + '/mean_VD.npy'),
               'mean_AV_GFP': np.load(res_folder + '/mean_AV.npy'),
               'mean_R_GFP': np.load(res_folder + '/mean_R.npy'),
               'distances_AD': np.load(res_folder + '/distances_AD.npy'),
               'distances_VD': np.load(res_folder + '/distances_VD.npy'),
               'distances_AV': np.load(res_folder + '/distances_AV.npy'),
               'distances_R': np.load(res_folder + '/distances_R.npy')}
        # also add sprouting front to the dictionary, if it was computed
        if os.path.exists(res_folder + '/distances_SF.npy'):
            #print('Sprouting front mask exists!!')
            res['distances_SF'] = np.load(res_folder + '/distances_SF.npy')

        res_list.append(res)

    # restore np.load for future normal usage
    #np.load = np_load_old
    np.load.__defaults__=(None, False, True, 'ASCII')

    df = pd.DataFrame(res_list)

    df['numberOfPixels_R'] = [len(df['distances_R'][i]) for i in range(len(df['distances_R']))]
    df['numberOfPixels_AV'] = [len(df['distances_AV'][i]) for i in range(len(df['distances_AV']))]
    
    if 'distances_SF' in df.keys():
        df['numberOfPixels_inSF'] = [len(np.where(df['distances_SF'][i]==0)[0]) for i in range(len(df['distances_SF']))]
        df['rationSF_to_totalNumber'] = [len(np.where(df['distances_SF'][i]==0)[0])/len(df['distances_AV'][i]) for i in range(len(df['distances_SF']))]
        
    
    df['TP'] = [TP_from_filename(fn_) for fn_ in df['retina_name']]
    df['ExperimentID'] = [expID_from_results_path(fn_) for fn_ in df['retina_name']]
    df['GFP'] = [1 if get_image_type_from_path(fn_)=='GFP' else 0 for fn_ in df['retina_name']]
    df['GFP X ERG'] = [1 if get_image_type_from_path(fn_)=='GFP X ERG' else 0 for fn_ in df['retina_name']]
    df['ERG'] = [1 if get_image_type_from_path(fn_)=='ERG' else 0 for fn_ in df['retina_name']]

    return df

def rearrange_results_for_plotting(parameters, df):
    # initiate struct for results
    #all_data = {'P5_P6': {}, 'P5_P7': {}, 'P5_P8': {}, 'P5_P9': {}, 'P5_P11': {}, 'P5_P12': {}, 'P8_P9': {}}
    all_data = {tp:{} for tp in df['TP'].unique()}

    for tp in all_data.keys():
        print('TIME PPOINT: ' + tp)
        all_data[tp]['color'] = assign_color(tp)
        all_data[tp]['injection_day'] = get_injectionDay(tp)
        all_data[tp]['culling_day'] = get_cullingDay(tp)
        all_data[tp]['migration_time'] = get_migrationTime(tp)
        
        all_data[tp]['mean_AD_GFP'] = np.zeros(np.shape(df['hist_AD_GFP'].iloc[0][0]))
        all_data[tp]['var_AD_GFP'] = np.zeros(np.shape(df['hist_AD_GFP'].iloc[0][0]))
        all_data[tp]['mean_AD_GFP_norm'] = np.zeros(np.shape(df['hist_AD_GFP'].iloc[0][0]))
        all_data[tp]['mean_VD_GFP'] = np.zeros(np.shape(df['hist_VD_GFP'].iloc[0][0]))
        all_data[tp]['var_VD_GFP'] = np.zeros(np.shape(df['hist_VD_GFP'].iloc[0][0]))
        all_data[tp]['mean_VD_GFP_norm'] = np.zeros(np.shape(df['hist_VD_GFP'].iloc[0][0]))
        all_data[tp]['mean_R_GFP'] = np.zeros(np.shape(df['hist_R_GFP'].iloc[0][0]))
        all_data[tp]['var_R_GFP'] = np.zeros(np.shape(df['hist_R_GFP'].iloc[0][0]))
        all_data[tp]['mean_R_GFP_norm'] = np.zeros(np.shape(df['hist_R_GFP'].iloc[0][0]))
        all_data[tp]['mean_AV_GFP'] = np.zeros(np.shape(df['hist_AV_GFP'].iloc[0][0]))
        all_data[tp]['var_AV_GFP'] = np.zeros(np.shape(df['hist_AV_GFP'].iloc[0][0]))
        all_data[tp]['mean_AV_R_GFP'] = np.zeros(np.shape(df['hist_AV_R_GFP'].iloc[0][0]))
        all_data[tp]['mean_AV_R_GFP_norm'] = np.zeros(np.shape(df['hist_AV_R_GFP'].iloc[0][0]))
        all_data[tp]['mean_AV_GFP_norm'] = np.zeros(np.shape(df['hist_AV_GFP'].iloc[0][0]))
        all_data[tp]['v_AV_mps'] = np.zeros(np.shape(df['v_AV_mps'].iloc[0]))
        all_data[tp]['v_R_mps'] = np.zeros(np.shape(df['v_R_mps'].iloc[0]))
        all_data[tp]['v_AD_mps'] = np.zeros(np.shape(df['v_AD_mps'].iloc[0]))
        all_data[tp]['v_VD_mps'] = np.zeros(np.shape(df['v_VD_mps'].iloc[0]))

        all_data[tp]['distances_AD'] = np.array([])
        all_data[tp]['distances_VD'] = np.array([])
        all_data[tp]['distances_AV'] = np.array([])
        all_data[tp]['distances_R'] = np.array([])
        
        if parameters['radial_bins'] != []:
            for i in parameters['radial_bins']:
                all_data[tp]['distances_AD_' + str(i[0]) + 'to' + str(i[1])] = np.array([])
                all_data[tp]['distances_VD_' + str(i[0]) + 'to' + str(i[1])] = np.array([])
                all_data[tp]['distances_AV_' + str(i[0]) + 'to' + str(i[1])] = np.array([])
                all_data[tp]['distances_R_' + str(i[0]) + 'to' + str(i[1])] = np.array([])

        all_data[tp]['mean_artery_cells'] = 0
        all_data[tp]['mean_vein_cells'] = 0
        all_data[tp]['mean_artery_vein_cell_ratio'] = 0
        
        if 'distances_SF' in df.keys():
            all_data[tp]['distances_AV_inSF'] = np.array([])
            all_data[tp]['distances_AV_outSF'] = np.array([])
            all_data[tp]['distances_AD_inSF'] = np.array([])
            all_data[tp]['distances_AD_outSF'] = np.array([])
            all_data[tp]['distances_VD_inSF'] = np.array([])
            all_data[tp]['distances_VD_outSF'] = np.array([])
            all_data[tp]['distances_R_inSF'] = np.array([])
            all_data[tp]['distances_R_outSF'] = np.array([])

        retina_count = 0

        for i, retina in df[df['TP'] == tp].iterrows():
            print('   Experiment ID: %s' %retina['ExperimentID'])
            retina_count += 1
            all_data[tp]['mean_AD_GFP'] += retina['hist_AD_GFP'][0]/np.sum(retina['hist_AD_GFP'][0])
            all_data[tp]['mean_VD_GFP'] += retina['hist_VD_GFP'][0]/np.sum(retina['hist_VD_GFP'][0])
            all_data[tp]['mean_AV_GFP'] += retina['hist_AV_GFP'][0]/np.sum(retina['hist_AV_GFP'][0])
            all_data[tp]['mean_AV_R_GFP'] += retina['hist_AV_R_GFP'][0]/np.sum(retina['hist_AV_R_GFP'][0])
            all_data[tp]['mean_R_GFP'] += retina['hist_R_GFP'][0]/np.sum(retina['hist_R_GFP'][0])
            all_data[tp]['v_AV_mps'] = retina['v_AV_mps']
            all_data[tp]['v_R_mps'] = retina['v_R_mps']
            all_data[tp]['v_AD_mps'] = retina['v_AD_mps']
            all_data[tp]['v_VD_mps'] = retina['v_VD_mps']

            if len(retina['distances_AD']) < parameters['pixel_no_per_retina']:
                all_data[tp]['distances_AD'] = np.concatenate((all_data[tp]['distances_AD'], retina['distances_AD']))
                all_data[tp]['distances_VD'] = np.concatenate((all_data[tp]['distances_VD'], retina['distances_VD']))
                all_data[tp]['distances_AV'] = np.concatenate((all_data[tp]['distances_AV'], retina['distances_AV']))
                all_data[tp]['distances_R'] = np.concatenate((all_data[tp]['distances_R'], retina['distances_R']))
                print("   Warning! Retina had too few pixels to draw randomly. Desired: %s, Available: %s" %(parameters['pixel_no_per_retina'],len(retina['distances_AD'])))
            else:
                all_data[tp]['distances_AD'] = np.concatenate((all_data[tp]['distances_AD'], np.random.choice(retina['distances_AD'], size=parameters['pixel_no_per_retina'], replace=False)))
                all_data[tp]['distances_VD'] = np.concatenate((all_data[tp]['distances_VD'], np.random.choice(retina['distances_VD'], size=parameters['pixel_no_per_retina'], replace=False)))
                all_data[tp]['distances_AV'] = np.concatenate((all_data[tp]['distances_AV'], np.random.choice(retina['distances_AV'], size=parameters['pixel_no_per_retina'], replace=False)))
                all_data[tp]['distances_R'] = np.concatenate((all_data[tp]['distances_R'], np.random.choice(retina['distances_R'], size=parameters['pixel_no_per_retina'], replace=False)))

            all_data[tp]['mean_vein_cells'] += np.sum(retina['hist_AV_GFP'][0][0:5])
            all_data[tp]['mean_artery_cells'] += np.sum(retina['hist_AV_GFP'][0][-5:])
            all_data[tp]['mean_artery_vein_cell_ratio'] += all_data[tp]['mean_vein_cells'] / all_data[tp]['mean_artery_cells']
            
            # only take values within a certain radial distance range
            if parameters['radial_bins'] != []:
                print('      Organizing data for different radial bins...')
                for l in parameters['radial_bins']:
                    indices = get_radial_indices(retina['distances_R'], l[0], l[1])
                    try:
                        indices = np.random.choice(indices, size=parameters['pixel_no_per_retina'], replace=False)
                    except:
                        print('      Warning! Radial bin %s: Too few pixels to draw randomly. Desired: %s, Available: %s' %(l,parameters['pixel_no_per_retina'],len(indices)))
                    all_data[tp]['distances_AV_' + str(l[0]) + 'to' + str(l[1])] = np.concatenate((all_data[tp]['distances_AV_' + str(l[0]) + 'to' + str(l[1])],[retina['distances_AV'][i] for i in indices]))
                    all_data[tp]['distances_AD_' + str(l[0]) + 'to' + str(l[1])] = np.concatenate((all_data[tp]['distances_AD_' + str(l[0]) + 'to' + str(l[1])],[retina['distances_AD'][i] for i in indices]))
                    all_data[tp]['distances_VD_' + str(l[0]) + 'to' + str(l[1])] = np.concatenate((all_data[tp]['distances_VD_' + str(l[0]) + 'to' + str(l[1])],[retina['distances_VD'][i] for i in indices]))
                    all_data[tp]['distances_R_' + str(l[0]) + 'to' + str(l[1])] = np.concatenate((all_data[tp]['distances_R_' + str(l[0]) + 'to' + str(l[1])], [retina['distances_R'][i] for i in indices]))
            
            if 'distances_SF' in df.keys():
                print('      Splitting up data according to sprouting front...')
                indices_inSF = get_radial_indices(retina['distances_SF'], 0, 0)
                try:
                    indices_inSF = np.random.choice(indices_inSF, size=parameters['pixel_no_per_retina'], replace=False)
                except:
                    print('      Warning! Too few pixels to draw randomly (inside SF). Desired: %s, Available: %s' %(parameters['pixel_no_per_retina'],len(indices_inSF)))
                
                all_data[tp]['distances_AV_inSF'] = np.concatenate((all_data[tp]['distances_AV_inSF'],[retina['distances_AV'][i] for i in indices_inSF]))
                all_data[tp]['distances_AD_inSF'] = np.concatenate((all_data[tp]['distances_AD_inSF'],[retina['distances_AD'][i] for i in indices_inSF]))
                all_data[tp]['distances_VD_inSF'] = np.concatenate((all_data[tp]['distances_VD_inSF'],[retina['distances_VD'][i] for i in indices_inSF]))
                all_data[tp]['distances_R_inSF'] = np.concatenate((all_data[tp]['distances_R_inSF'],[retina['distances_R'][i] for i in indices_inSF]))
                
                indices_outSF = get_radial_indices(retina['distances_SF'], 150, 1500)
                try:
                    indices_outSF = np.random.choice(indices_outSF, size=parameters['pixel_no_per_retina'], replace=False)
                except:
                    print('      Warning! Too few pixels to draw randomly (outside SF). Desired: %s, Available: %s' %(parameters['pixel_no_per_retina'],len(indices_outSF)))
            
                all_data[tp]['distances_AV_outSF'] = np.concatenate((all_data[tp]['distances_AV_outSF'],[retina['distances_AV'][i] for i in indices_outSF]))
                all_data[tp]['distances_AD_outSF'] = np.concatenate((all_data[tp]['distances_AD_outSF'],[retina['distances_AD'][i] for i in indices_outSF]))
                all_data[tp]['distances_VD_outSF'] = np.concatenate((all_data[tp]['distances_VD_outSF'],[retina['distances_VD'][i] for i in indices_outSF]))
                all_data[tp]['distances_R_outSF'] = np.concatenate((all_data[tp]['distances_R_outSF'],[retina['distances_R'][i] for i in indices_outSF]))
            
        all_data[tp]['retina_count'] = retina_count

        all_data[tp]['mean_AD_GFP'] /= retina_count
        all_data[tp]['mean_VD_GFP'] /= retina_count
        all_data[tp]['mean_R_GFP'] /= retina_count
        all_data[tp]['mean_AV_GFP'] /= retina_count
        all_data[tp]['mean_AV_R_GFP'] /= retina_count
        all_data[tp]['mean_artery_cells'] /= retina_count
        all_data[tp]['mean_vein_cells'] /= retina_count
        all_data[tp]['mean_artery_vein_cell_ratio'] /= retina_count
        
        all_data[tp]['median_AV'] = np.median(all_data[tp]['distances_AV'])
        all_data[tp]['median_AD'] = np.median(all_data[tp]['distances_AD'])
        all_data[tp]['median_VD'] = np.median(all_data[tp]['distances_VD'])
        all_data[tp]['median_R'] = np.median(all_data[tp]['distances_R'])
        all_data[tp]['mean_AV'] = np.mean(all_data[tp]['distances_AV'])
        all_data[tp]['mean_R'] = np.mean(all_data[tp]['distances_R'])
        all_data[tp]['mean_AD'] = np.mean(all_data[tp]['distances_AD'])
        all_data[tp]['mean_VD'] = np.mean(all_data[tp]['distances_VD'])
        all_data[tp]['std_AV'] = np.std(all_data[tp]['distances_AV'])
        all_data[tp]['std_R'] = np.std(all_data[tp]['distances_R'])
        all_data[tp]['std_AD'] = np.std(all_data[tp]['distances_AD'])
        all_data[tp]['std_VD'] = np.std(all_data[tp]['distances_VD'])
        
        for i in parameters['radial_bins']:
            all_data[tp]['median_AV_' + str(i[0]) + 'to' + str(i[1])] = np.median(all_data[tp]['distances_AV_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['median_AD_' + str(i[0]) + 'to' + str(i[1])] = np.median(all_data[tp]['distances_AD_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['median_VD_' + str(i[0]) + 'to' + str(i[1])] = np.median(all_data[tp]['distances_VD_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['median_R_' + str(i[0]) + 'to' + str(i[1])] = np.median(all_data[tp]['distances_R_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['mean_AV_' + str(i[0]) + 'to' + str(i[1])] = np.mean(all_data[tp]['distances_AV_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['mean_R_' + str(i[0]) + 'to' + str(i[1])] = np.mean(all_data[tp]['distances_R_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['mean_AD_' + str(i[0]) + 'to' + str(i[1])] = np.mean(all_data[tp]['distances_AD_' + str(i[0]) + 'to' + str(i[1])])
            all_data[tp]['mean_VD_' + str(i[0]) + 'to' + str(i[1])] = np.mean(all_data[tp]['distances_VD_' + str(i[0]) + 'to' + str(i[1])])
        
        if 'distances_SF' in df.keys():
            all_data[tp]['median_AV_inSF'] = np.median(all_data[tp]['distances_AV_inSF'])
            all_data[tp]['median_AD_inSF'] = np.median(all_data[tp]['distances_AD_inSF'])
            all_data[tp]['median_VD_inSF'] = np.median(all_data[tp]['distances_VD_inSF'])
            all_data[tp]['median_AV_outSF'] = np.median(all_data[tp]['distances_AV_outSF'])
            all_data[tp]['median_AD_outSF'] = np.median(all_data[tp]['distances_AD_outSF'])
            all_data[tp]['median_VD_outSF'] = np.median(all_data[tp]['distances_VD_outSF'])

    return pd.DataFrame(all_data)

def confidence_KS_test(parameters, df, tp, exp_conditions):
    """
    Do the KS test a specified number of times to get a confidence interval of the p-value.
    """
    typesOfDistances = ['AD', 'VD', 'AV', 'R']
    """
    
    # prepare the dictionary 
    all_data = {channel:{} for channel in parameters['image_types']}
    p_values = {channel:{} for channel in parameters['image_types']}
    
    for channel in parameters['image_types']:
        all_data[channel] = {cond:{} for cond in parameters['load_conditions']}
        
        for cond in parameters['load_conditions']:
            all_data[channel][cond] = {tp:{} for tp in df['TP'].unique()}
    
        for number in range(parameters['number_ks_test']):
            # initialize everything we need

            for typeOfDistance in typesOfDistances:
                all_data[channel][cond][tp]['distances_' + typeOfDistance] = np.array([])
            
                if parameters['radial_bins'] != []:
                    for i in parameters['radial_bins']:
                        all_data[channel][cond][tp]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])] = np.array([])
                
                if 'distances_SF' in df.keys():              
                    all_data[channel][cond][tp]['distances_' + typeOfDistance + '_inSF'] = np.array([])
                    all_data[channel][cond][tp]['distances_' + typeOfDistance + '_outSF'] = np.array([])
    """
    p_values = {channel:{} for channel in parameters['image_types']}
    rearranged_data = {}
    for image_type in parameters['image_types']:
        p_values[image_type] = {}
            
        for typeOfDistance in typesOfDistances:
            p_values[image_type]['distances_' + typeOfDistance] = np.array([])
            
            if parameters['radial_bins'] != []:
                for i in parameters['radial_bins']:
                    p_values[image_type]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])] = np.array([])
                
            if 'distances_SF' in df[image_type].keys():              
                p_values[image_type]['distances_' + typeOfDistance + '_inSF'] = np.array([])
                p_values[image_type]['distances_' + typeOfDistance + '_outSF'] = np.array([])
        
        for number in range(parameters['number_ks_test']):
            print('')
            print('KS TEST ROUND ' + str(number))
            print('')
            rearranged_data[image_type] = {}
            for i in exp_conditions:
                rearranged_data[image_type][i] = rearrange_results_for_plotting(parameters,df[image_type][df[image_type][i]==1])
            
            for typeOfDistance in typesOfDistances:
                p_values[image_type]['distances_' + typeOfDistance] = np.append(p_values[image_type]['distances_' + typeOfDistance], scipy.stats.kstest(rearranged_data[image_type][exp_conditions[0]][tp]['distances_' + typeOfDistance], rearranged_data[image_type][exp_conditions[1]][tp]['distances_' + typeOfDistance], mode='exact')[1])
            
                if parameters['radial_bins'] != []:
                    for i in parameters['radial_bins']:
                        p_values[image_type]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])] = np.append(p_values[image_type]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])], scipy.stats.kstest(rearranged_data[image_type][exp_conditions[0]][tp]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])], rearranged_data[image_type][exp_conditions[1]][tp]['distances_' + typeOfDistance + '_' + str(i[0]) + 'to' + str(i[1])], mode='exact')[1])

                if 'distances_SF' in df[image_type].keys():    
                    p_values[image_type]['distances_' + typeOfDistance + '_inSF'] = np.append(p_values[image_type]['distances_' + typeOfDistance + '_inSF'], scipy.stats.kstest(rearranged_data[image_type][exp_conditions[0]][tp]['distances_' + typeOfDistance + '_inSF'], rearranged_data[image_type][exp_conditions[1]][tp]['distances_' + typeOfDistance + '_inSF'], mode='exact')[1])
                    p_values[image_type]['distances_' + typeOfDistance + '_outSF'] = np.append(p_values[image_type]['distances_' + typeOfDistance + '_outSF'], scipy.stats.kstest(rearranged_data[image_type][exp_conditions[0]][tp]['distances_' + typeOfDistance + '_outSF'], rearranged_data[image_type][exp_conditions[1]][tp]['distances_' + typeOfDistance + '_outSF'], mode='exact')[1])

    return p_values

def job_finished():
    print()
    print('DOOONEE!!!') 
    
    lucky_number = random.randint(0, 6)
    
    if lucky_number == 0:
        print('Also, as we are too busy to make it to the Naturkundemuseum - here is a dinosaur:')
        print('                                              ____')
        print('   ___                                      .-~. /_"-._')
        print('  `-._~-.                                  / /_ "~o\  :o')
        print("      \  \                                / : \~x.  ` ')")
        print('       ]  Y                              /  |  Y< ~-.__j')
        print('      /   !                        _.--~T : l  l<  /.-~')
        print('     /   /                 ____.--~ .   ` l /~\ \<|Y' )
        print('    /   /             .-~~"' + "        /| .    ',-~\ \L|")
        print('   /   /             /     .^   \ Y~Y \.^>/l_   "' + "--'")
        print('  /   Y           .-"(  .  l__  j_j l_/ /~_.-~    .')
        print(' Y    l          /    \  )    ~~~." / `/"~ / \.__/l_')
        print(' |     \     _.-"      ~-{__     l  :  l._Z~-.___.--~')
        print(' |      ~---~           /   ~~"---\_  ' + "' __[>")
        print(' l  .                _.^   ___     _>-y~')
        print('  \  \     .      .-~   .-~   ~>--"  /')
        print('   \  ~---"            /     ./  _.-' + "'")
        print('    "-.,_____.,_  _.--~\     _.-~')
        print('                ~~     (   _}')
        print("                        `. ~(")
        print("                          )  \ ")
        print("                         /,`--'~\--' ")
    
    if lucky_number == 1:
        print("On top, here is a random anteater. You're welcome.")
        print('       _.---._    /\\')
        print("    ./'       " + '"--`\//   ')
        print('  ./              o \          .-----.')
        print(' /./\  )______   \__ \        ( help! )')
        print("./  / /\ \   | \ \  \ \       /`-----'")
        print('   / /  \ \  | |\ \  \7--- ooo ooo ooo ooo ooo ooo')

    if lucky_number == 2:    
        print('Now, what about some impossible forks?')
        print(' __________________         __________________')
        print('()_________________\       /_________________()')
        print(' ________________  |       |  ________________')
        print('()____________| |  |       |  | |____________()')
        print(' ______________\|  |       |  |/______________')
        print('()_________________|       |_________________()')
        print()
        print('WOOOOOOOOOW')
    
    if lucky_number == 3:
        print("Also, here's a maze for you. Cheers!")
        print(" __________________________________   ")
        print("| _____ |   | ___ | ___ ___ | |   | |")
        print('| |   | |_| |__ | |_| __|____ | | | |')
        print('| | | |_________|__ |______ |___|_| |')
        print('| |_|   | _______ |______ |   | ____|')
        print('| ___ | |____ | |______ | |_| |____ |')
        print('|___|_|____ | |   ___ | |________ | |')
        print('|  _________| | |__ | |______ | | | |')
        print('| | | ________| | __|____ | | | __| |')
        print('|_| |__ |   | __|__ | ____| | |_| __|')
        print('|  _____| | |_____| |__|    |__ |__ |')
        print('| |_______|_______|___|___|___|_____|')

    if lucky_number == 4:
        print('                    .+~                :xx++::')
        print('                   :`. -          .!!X!~"?!`~!~!. :-:.')
        print('                  {             .!!!H":.~ ::+!~~!!!~ `%X.')
        print("                  '             ~~!M!!>!!X?!!!!!!!!!!...!~.")
        print('                              {!:!MM!~:XM!!!!!!.:!..~ !.  `{')
        print('                  {: `   :~ .:{~!!M!XXHM!!!X!XXHtMMHHHX!  ~ ~')
        print("                ~~~~{' ~!!!:!!!!!XM!!M!!!XHMMMRMSXXX!!!!!!:  {`")
        print('                  `{  {::!!!!!X!X?M!!M!!XMMMMXXMMMM??!!!!!?!:~{')
        print("               : '~~~{!!!XMMH!!XMXMXHHXXXXM!!!!MMMMSXXXX!!!!!!!~")
        print('            :    ::`~!!!MMMMXXXtMMMMMMMMMMMHX!!!!!!HMMMMMX!!!!!: ~')
        print("               '~:~!!!!!MMMMMMMMMMMMMMMMMMMMMMXXX!!!M??MMMM!!X!!i:")
        print('               {~{!!!!!XMMMMMMMMMMMM8M8MMMMM8MMMMMXX!!!!!!!!X!?t?!:')
        print('               ~:~~!!!!?MMMMMM@M@RMRRR$@@MMRMRMMMMMMXSX!!!XMMMX{?X!')
        print('             :XX {!!XHMMMM88MM88BR$M$$$$8@8RN88MMMMMMMMHXX?MMMMMX!!!')
        print('           .:X! {XMSM8M@@$$$$$$$$$$$$$$$$$$$B8R$8MMMMMMMMMMMMMMMMX!X')
        print('          :!?! !?XMMMMM8$$$$8$$$$$$$$$$$$$$BBR$$MMM@MMMMMMMMMMMMMM!!X')
        print('        ~{!!~ {!!XMMMB$$$$$$$$$$$$$$$$$$$$$$$$MMR$8MR$MMMMMMMMMMMMM!?!:')
        print('        :~~~ !:X!XMM8$$$$$$$$$$$$$$$$$$$$$$$RR$$MMMMR8NMMMMMMMMMMMMM{!`-')
        print("    ~:{!:~`~':!:HMM8N$$$$$$$$$$$$$$$$$$$$$$$$$8MRMM8R$MRMMMMMMMMRMMMX!")
        print('  !X!``~~   :~XM?SMM$B$$$$$$$$$$$$$$$$$$$$$$BR$$MMM$@R$M$MMMMMM$MMMMX?L')
        print(' X~.      : `!!!MM#$RR$$$$$$$$$$$$$$$$$R$$$$$R$M$MMRRRM8MMMMMMM$$MMMM!?:')
        print(' ! ~ {~  !! !!~`` :!!MR$$$$$$$$$$RMM!?!??RR?#R8$M$MMMRM$RMMMM8MM$MMM!M!:>')
        print(": ' >!~ '!!  !   .!XMM8$$$$$@$$$R888HMM!!XXHWX$8$RM$MR5$8MMMMR$$@MMM!!!{ ~")
        print("!  ' !  ~!! :!:XXHXMMMR$$$$$$$$$$$$$$$$8$$$$8$$$MMR$M$$$MMMMMM$$$MMM!!!!")
        print(' ~{!!!  !!! !!HMMMMMMMM$$$$$$$$$$$$$$$$$$$$$$$$$$MMM$M$$MM8MMMR$$MMXX!!!!/:`')
        print('  ~!!!  !!! !XMMMMMMMMMMR$$$$$$$$$$$$R$RRR$$$$$$$MMMM$RM$MM8MM$$$M8MMMX!!!!:')
        print('  !~ ~  !!~ XMMM%!!!XMMX?M$$$$$$$$B$MMSXXXH?MR$$8MMMM$$@$8$M$B$$$$B$MMMX!!!!')
        print("  ~!    !! 'XMM?~~!!!MMMX!M$$$$$$MRMMM?!%MMMH!R$MMMMMM$$$MM$8$$$$$$MR@M!!!!!")
        print('  {>    !!  !Mf x@#"~!t?M~!$$$$$RMMM?Xb@!~`??MS$M@MMM@RMRMMM$$$$$$RMMMMM!!!!')
        print("  !    '!~ {!!:!?M   !@!M{XM$$R5M$8MMM$! -XXXMMRMBMMM$RMMM@$R$BR$MMMMX??!X!!")
        print("  !    '!  !!X!!!?::xH!HM:MM$RM8M$RHMMMX...XMMMMM$RMMRRMMMMMMM8MMMMMMMMX!!X!")
        print('  !     ~  !!?:::!!!MXMR~!MMMRMM8MMMMMS!!M?XXMMMMM$$M$M$RMMMM8$RMMMMMMMM%X!!')
        print('  ~     ~  !~~X!!XHMMM?~ XM$MMMMRMMMMMM@MMMMMMMMMM$8@MMMMMMMMRMMMMM?!MMM%HX!')
        print('           !!!!XSMMXXMM .MMMMMMMM$$$BB8MMM@MMMMMMMR$RMMMMMMMMMMMMMMMXX!?H!XX')
        print('           XHXMMMMMMMM!.XMMMMMMMMMR$$$8M$$$$$M@88MMMMMMMMMMMMMMM!XMMMXX!!!XM')
        print('      ~   {!MMMMMMMMRM:XMMMMMMMMMM8R$$$$$$$$$$$$$$$NMMMMMMMM?!MM!M8MXX!!/t!M')
        print("      '   ~HMMMMMMMMM~!MM8@8MMM!MM$$8$$$$$$$$$$$$$$8MMMMMMM!!XMMMM$8MR!MX!MM")
        print("          'MMMMMMMMMM'MM$$$$$MMXMXM$$$$$$$$$$$$$$$$RMMMMMMM!!MMM$$$$MMMMM{!M")
        print("          'MMMMMMMMM!'MM$$$$$RMMMMMM$$$$$$$$$$$$$$$MMM!MMMX!!MM$$$$$M$$M$M!M")
        print('           !MMMMMM$M! !MR$$$RMM8$8MXM8$$$$$$$$$$$$NMMM!MMM!!!?MRR$$RXM$$MR!M')
        print('           !M?XMM$$M.{ !MMMMMMSUSRMXM$8R$$$$$$$$$$#$MM!MMM!X!t8$M$MMMHMRMMX$')
        print("    ,-,   '!!!MM$RMSMX:.?!XMHRR$RM88$$$8M$$$$$R$$$$8MM!MMXMH!M$$RMMMMRNMMX!$")
        print("   -'`    '!!!MMMMMMMMMM8$RMM8MBMRRMR8RMMM$$$$8$8$$$MMXMMMMM!MR$MM!M?MMMMMM$")
        print("          'XX!MMMMMMM@RMM$MM@$$BM$$$M8MMMMR$$$$@$$$$MM!MMMMXX$MRM!XH!!??XMMM")
        print('          `!!!M?MHMMM$RMMMR@$$$$MR@MMMM8MMMM$$$$$$$WMM!MMMM!M$RMM!!.MM!%M?~!')
        print('           !!!!!!MMMMBMM$$RRMMMR8MMMMMRMMMMM8$$$$$$$MM?MMMM!f#RM~    `~!!!~!')
        print('           ~!!HX!!~!?MM?MMM??MM?MMMMMMMMMRMMMM$$$$$MMM!MMMM!!')
        print("           '!!!MX!:`~~`~~!~~!!!!XM!!!?!?MMMM8$$$$$MMMMXMMM!!")
        print('            !!~M@MX.. {!!X!!!!XHMHX!!``!XMMMB$MM$$B$M!MMM!!')
        print('            !!!?MRMM!:!XHMHMMMMMMMM!  X!SMMX$$MM$$$RMXMMM~')
        print('             !M!MMMM>!XMMMMMMMMXMM!!:!MM$MMMBRM$$$$8MMMM~')
        print("             `?H!M$R>'MMMM?MMM!MM6!X!XM$$$MM$MM$$$$MX$f")
        print('              `MXM$8X MMMMMMM!!MM!!!!XM$$$MM$MM$$$RX@"')
        print('               ~M?$MM !MMMMXM!!MM!!!XMMM$$$8$XM$$RM!`')
        print('                !XMMM !MMMMXX!XM!!!HMMMM$$$$RH$$M!~')
        print("                'M?MM `?MMXMM!XM!XMMMMM$$$$$RM$$#")
        print('                 `>MMk ~MMHM!XM!XMMM$$$$$$BRM$M"')
        print('                  ~`?M. !M?MXM!X$$@M$$$$$$RMM#')
        print('                    `!M  !!MM!X8$$$RM$$$$MM#`')
        print('                      !% `~~~X8$$$$8M$$RR#`')
        print('                       !!x:xH$$$$$$$R$R*`')
        print('                        ~!?MMMMRRRM@M#`     ') 
        print('                         `~???MMM?M"`')
        print('                             ``~~')
    
    if lucky_number == 5:                             
        print('                 /`.    /`.')
        print("                f   \  ,f  \ ")
        print('    Gee Brain,  |    \/-`\  \      The same thing we do')
        print(' what do you    i.  _\';.,X j      every night, Pinky.')
        print('   want to do    `:_\ (  \ \',-.   Try to take over')
        print("        tonight?   .'" + '"`\ a\eY' +"' )   the world!  _,. ")
        print('                   `._"\`-' + "' `-/            .-;'  |")
        print("                     /;-`._.-';\.        ,'" + ',"    |')
        print("                   .'/   "+'     | `\.-'"" + '""-/ /      j')
        print('                 ,/ /         i,-"        (  ,/  /')
        print("              .-' .f         .'            `" + '"/  /')
        print('             / ,,/ffj\      /          .-"`.' + "-.'")
        print("            / /_\`--//)     \ ,--._ .-'_,-'; /")
        print('           f  ".-"-._;' + "'      `._ _.,-i; /_; /")
        print("           `.,'   |; \          \`\_,/-'  \'")
        print("            .'    l \ `.        /" + '"\ _ \`  j')
        print("            f      : `-'        `._;." + '"/`-' + "'")
        print('            |      `.               ,7  \ ')
        print("            l       j             .'/ - \`.")
        print("           .j.  .   <            (.'    .\ \f`. |\,'")
        print("          ,' `.  \ / \           `|      \,'||-:j")
        print("        .'  .'\   Y.  \___......__\ ._   /`.||")
        print('__.._,-" .-"' + "'" + '"")  /' ,' _          \ |  /"-.`j""``---.._')
        print("  .'_.-'" + '"     / .("-' + "'-" + '"":\        ._)|_(__. "')
        print(" ;.'         /-'---" + '"".--"' +"'      /,_,^-._ .)")
        print(" `:=.__.,itz `---._.;' ")
        
    if lucky_number == 6:
        print("Now it's time to relax :)")
        print('                   \       /            _\/_')
        print("                     .-'-.              //o\  _\/_")
        print('  _  ___  __  _ --_ /     \ _--_ __  __ _ | __/o\\ _')
        print("=-=-_=-=-_=-=_=-_= -=======- = =-=_=-=_,-'|''''-|-,_ ")
        print(' =- _=-=-_=- _=-= _--=====- _=-=_-_,-"          |')
        print('jgs=- =- =-= =- = -  -===- -= - ."        ')
