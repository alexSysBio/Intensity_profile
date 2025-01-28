# -*- coding: utf-8 -*-
"""

Alexandros Papagiannakis
HHMI at Stanford University, Christine Jacobs-Wagner lab, Sarafan ChEM-H, 2025

"""

import nd2_to_array as ndtwo
import Biviriate_medial_axis_estimation as medax
from skimage.filters import threshold_otsu
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops_table
import pandas as pd
import numpy as np
import os

def load_images(images_path):
    return ndtwo.nd2_to_array(images_path)


def get_specific_image(ndtwo_object, xy_position, channel, timepoint):
    return ndtwo_object[2][xy_position-1][channel][timepoint]


def get_otsu_mask(phase_image, min_area, max_area):
    masked_image = phase_image >  threshold_otsu(phase_image)
    masked_image = ~ masked_image
    
    labeled_image = label(masked_image)
    properties = ['label', 'area']
    region_table = regionprops_table(labeled_image, properties=properties)
    region_df= pd.DataFrame(region_table)
    region_df = region_df[region_df.area.between(min_area, max_area)]
    
    return label(np.isin(labeled_image, region_df.label))


def plot_cell_labels(labeled_image, save_suffix, save_path):
    
    properties = ['label', 'area', 'centroid']
    label_df = pd.DataFrame(regionprops_table(labeled_image, properties=properties))
    
    plt.figure(figsize=(10,10))
    plt.title(save_suffix, fontsize=14, style='italic')
    plt.imshow(labeled_image)
    plt.xticks([])
    plt.yticks([])
    for index, row in label_df.iterrows():
        cen_x = row['centroid-1']
        cen_y = row['centroid-0']
        lbl = row.label
        plt.text(cen_x, cen_y, f"{lbl:.0f}", color='white', fontsize=24)
    if os.path.isdir(save_path):
        plt.savefig(save_path+'/'+save_suffix+'_cell_labels.jpeg')
    plt.show()

def get_medial_axis_for_label(labeled_image, lbl, radius_px=8, half_angle=25, 
                                        cap_knot=16, max_degree=30, verbose=True):
    
    masked_image = labeled_image==lbl
    
    y_coords, x_coords = np.nonzero(masked_image)
    min_y = np.min(y_coords)-3
    max_y = np.max(y_coords)+4
    min_x = np.min(x_coords)-3
    max_x = np.max(x_coords)+4
    
    crop_pad = (min_x, min_y, max_x, max_y)
    
    cropped_cell_mask = masked_image[min_y:max_y, min_x:max_x]
    
    medial_axis_def = medax.get_medial_axis(cropped_cell_mask, radius_px, half_angle, 
                                            cap_knot, max_degree, verbose)
    oned_df = medax.get_oned_coordinates(cropped_cell_mask, medial_axis_def[0])
    
    return oned_df, crop_pad, cropped_cell_mask


def get_oned_intensity(crop_pad, oned_df, cropped_mask, fluor_image, channel):
    
    min_x, min_y, max_x, max_y = crop_pad
    cropped_fluor_image = fluor_image[min_y:max_y, min_x:max_x]
    oned_df['fluor_'+channel] = cropped_fluor_image[np.nonzero(cropped_mask)]
    
    return oned_df
    
    
def get_intensity_profiles(fluor_oned_df, width_zone=6, bin_number="auto"):
        
    if width_zone % 2 != 0:
        raise ValueError('The width_zone parameter should be an even integer')
    
    if bin_number == "auto":
        length_bins = np.arange(fluor_oned_df.arch_length.min(), fluor_oned_df.arch_length.max()+1, 1)
        fluor_oned_df['length_bin']= pd.cut(fluor_oned_df.arch_length, bins=length_bins)
    else:
        if type(bin_number) == int:
            fluor_oned_df['length_bin']= pd.cut(fluor_oned_df.arch_length, bins=bin_number)
        else:
            raise ValueError('The bin_number must be set to "auto" or be an integer')
    
    press_df = fluor_oned_df[fluor_oned_df.width.abs()<width_zone/2]
    mean_df = press_df.groupby('length_bin').mean()
    
    return mean_df
            
    
    
    



        
    