# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 14:31:47 2025

Alexandros Papagiannakis
HHMI at Stanford University, Christine Jacobs-Wagner lab, Sarafan ChEM-H, 2025

"""


import os
# sys.path.append(os.path.join(os.path.abspath(os.getcwd()), "/make_intensity_profile.py"))
import make_intensity_profile as mint
import matplotlib.pyplot as plt
from skimage.measure import find_contours
import uneven_background_correction as bkg
import numpy as np

# save_path= "/..._intensity_profiles"
# images_path = "/...nd2"
# xy_position = 3
# timepoint = 0
# images = mint.load_images(images_path)
# channels=['GFP', 'DAPI']
# channel_labels=['RplA-GFP', 'HupA-mCherry']                                                                                              



def get_save_suffix(xy_position, timepoint):
    return 'xy'+str(xy_position)+'_time'+str(timepoint)



def get_segmented_labels_and_images(images, xy_position, timepoint, save_path):
    
    channel_list = images[3]
    fluor_images_dict = {}

    for ch in channel_list:
        if ch == 'Phase' or ch == 'Trans':
            phase_image = mint.get_specific_image(images, xy_position, ch, timepoint)
        else:
            fluor_images_dict[ch] = mint.get_specific_image(images, xy_position, ch, timepoint)
    
    phase_labels = mint.get_otsu_mask(phase_image, min_area=20, max_area=5000)
    
    for ch in fluor_images_dict:
        bkg_cor_image = bkg.back_sub(fluor_images_dict[ch], phase_labels>0, 
                                     dilation=35, estimation_step=128, smoothing_sigma=60, show=False)
        fluor_images_dict[ch] = bkg_cor_image[0]
        
    save_suffix = get_save_suffix(xy_position, timepoint)
    mint.plot_cell_labels(phase_labels, save_suffix, save_path)
    
    return phase_labels, phase_image, fluor_images_dict
    


def select_label():
    i = 0
    while i == 0:
        try:
            lbl = int(input('Choose the label number of the desired cell:'))
        except ValueError:
            print('please choose a number')
            lbl = int(input('Choose the label number of the desired cell:'))
        i+=1
    return lbl



def check_if_crop_pad_is_within_bounds(crop_pad, square_sensor_dimension):
    min_x, min_y, max_x, max_y = crop_pad
    return min_x > 0 and min_y >0 and max_x < square_sensor_dimension and max_y < square_sensor_dimension


def get_mean_intensity_dataframe(cell_label, phase_labels, 
                                 bkg_cor_images_dict, save_path, medial_axis_threshold=20):
    
    if cell_label == 0:
        lbl = select_label()
    else:
        lbl = cell_label
        
    oned_df, crop_pad, cropped_cell_mask = mint.get_medial_axis_for_label(phase_labels, lbl, radius_px=8, half_angle=25, 
                                                                     cap_knot=16, max_degree=30, verbose=True)
    
    image_sensor = phase_labels.shape[0]
    
    cell_in_bounds = check_if_crop_pad_is_within_bounds(crop_pad, image_sensor)
    
    if oned_df.width.abs().max() <= medial_axis_threshold and cell_in_bounds == True:
    
        for ch in bkg_cor_images_dict:
            oned_df = mint.get_oned_intensity(crop_pad, oned_df, cropped_cell_mask, bkg_cor_images_dict[ch], ch)
        
        mean_int_df = mint.get_intensity_profiles(oned_df, width_zone=6, bin_number="auto")
        
        return mean_int_df, cropped_cell_mask, crop_pad, lbl
    
    elif oned_df.width.abs().max() > medial_axis_threshold:
        raise TypeError(f'This medial axis likely does not correspond to a real cell because it has a width above the specified width threshold of {medial_axis_threshold} pixels')
    elif cell_in_bounds == False:
        raise TypeError('This cell extends outside the field of view')


def plot_intensity_profile(mean_intensity_df, channels, channel_labels, save_suffix, save_path):
    
    
    fig, ax1 = plt.subplots()
    plt.title(save_suffix,  fontsize=14, style='italic')
    for i in range(len(channels)):
        if i == 0:
            ax1.plot(mean_intensity_df.scaled_length/2+0.5, mean_intensity_df['fluor_'+str(channels[i])], color='royalblue')
            ax1.set_ylabel(channel_labels[i], fontsize=14, fontweight='bold', color='royalblue')
        else:
            ax2 = ax1.twinx()
            ax2.plot(mean_intensity_df.scaled_length/2+0.5, mean_intensity_df['fluor_'+str(channels[i])], color='tomato')
            ax2.set_ylabel(channel_labels[i], fontsize=14, fontweight='bold', color='tomato')
    ax1.set_xlabel('Relative cell length', fontsize=14, fontweight='bold')
    
    if os.path.isdir(save_path):
        plt.savefig(save_path+'/intensity_profile_'+save_suffix+'.eps')
    plt.show()
    


def get_intensity_profile(phase_labels, bkg_cor_images_dict, channels, channel_labels, save_suffix, save_path):
    
    cell_label = 0
    mean_int_df, cropped_cell_mask, crop_pad, lbl = get_mean_intensity_dataframe(cell_label, phase_labels, 
                                                                                 bkg_cor_images_dict, save_path)
    
    labeled_save_suffix = save_suffix+'_label'+str(lbl)
    plot_intensity_profile(mean_int_df, channels, channel_labels, labeled_save_suffix, save_path)
    
    return mean_int_df, cropped_cell_mask, crop_pad, lbl, labeled_save_suffix



def plot_cell_images(phase_image, fluor_images_dict, cropped_cell_mask, crop_pad, save_suffix, save_path):
    
    cell_mesh_y, cell_mesh_x = zip(*find_contours(cropped_cell_mask, level=0.5)[0])
    
    min_x, min_y, max_x, max_y = crop_pad
    cropped_phase_image = phase_image[min_y:max_y, min_x:max_x]
    
    plt.imshow(cropped_phase_image, cmap='gray')
    plt.plot(cell_mesh_x, cell_mesh_y, color='yellow', linewidth=2)
    plt.plot([5,5+1/0.066], [5,5], linewidth=3, color='white')
    plt.text(5,20, '1um', fontsize=12, fontweight='bold', color='white')
    plt.xticks([])
    plt.yticks([])
    cbar = plt.colorbar()
    cbar.set_label('Phase contrast (a.u.)', rotation=90)
    if os.path.isdir(save_path):
        plt.savefig(save_path+'/'+save_suffix+'_phase_contrast_image.eps')
    plt.show()
    
    
    for ch in fluor_images_dict:
        cropped_fluor_image = fluor_images_dict[ch][min_y:max_y, min_x:max_x]
        plt.imshow(cropped_fluor_image, cmap='gray')
        plt.plot(cell_mesh_x, cell_mesh_y, color='yellow', linewidth=2)
        plt.plot([5,5+1/0.066], [5,5], linewidth=3, color='white')
        plt.text(5,20, '1um', fontsize=12, fontweight='bold', color='white')
        plt.xticks([])
        plt.yticks([])
        cbar = plt.colorbar()
        cbar.set_label(ch, rotation=90)
        if os.path.isdir(save_path):
            plt.savefig(save_path+'/'+save_suffix+'_'+ch+'_image.eps')
        plt.show()
    


def plot_intensity_profiles_and_cells(images, xy_position, timepoint, channels, channel_labels, save_path):
    
    phase_labels, phase_image, bkg_cor_images_dict = get_segmented_labels_and_images(images, xy_position, timepoint, save_path)
    save_suffix = get_save_suffix(xy_position, timepoint)
    
    mean_int_df, cropped_cell_mask, crop_pad, lbl, labeled_save_suffix = get_intensity_profile(phase_labels, bkg_cor_images_dict, 
                                                                                               channels, channel_labels, save_suffix, save_path)
    plot_cell_images(phase_image, bkg_cor_images_dict, cropped_cell_mask, crop_pad, labeled_save_suffix, save_path)



def get_intensity_profile_stats(images, xy_position, timepoint, save_path):
    
    phase_images_dict = {}
    channel_images_dict = {}
    mean_int_dict = {}
    
    channel_list = images[3]
    for ch in channel_list:
        if ch == 'Phase' or ch == 'Trans':
            phase_image = mint.get_specific_image(images, xy_position, ch, timepoint)
    
    phase_labels = mint.get_otsu_mask(phase_image, min_area=20, max_area=5000)
    cell_label_list = list(np.unique(phase_labels.ravel()))
    cell_label_list.remove(0)
    
    for cell_label in cell_label_list:
        try:
            mean_int_df, phase_image, fluor_images_dict, cropped_cell_mask, crop_pad, lbl = get_mean_intensity_dataframe(images, xy_position, timepoint, cell_label, save_path)
            phase_images_dict[cell_label] = phase_image
            channel_images_dict[cell_label] = fluor_images_dict
            mean_int_dict[cell_label] = mean_int_df
        except TypeError:
            print(f'Label {cell_label} is aborted because it does not correspond to a good segmentation instance...')
        except ValueError:
            print(f'Label {cell_label} out-of-bounds and aborted...')
    
    return mean_int_dict, phase_images_dict, channel_images_dict



# def get_intensity_profiles_for_many_frames(images, xy_positions = np.arange(1,10), time_points = [0,6,12,18,24,30,34,38,42]):
    

#     for xy_position in xy_positions:
#         for timepoint in time_points:
#             phase_labels = mint.get_otsu_mask(phase_image, min_area=20, max_area=5000)
#             get_intensity_profile_stats(images, xy_position, timepoint, cell_label_list, save_path)
            
            
#             for ch in channel_list:
#                 fluor_images_dict = {}
#                 if ch == 'Phase' or ch == 'Trans':
#                     phase_image = mint.get_specific_image(images, xy_position, ch, timepoint)
#                 else:
#                     fluor_images_dict[ch] = mint.get_specific_image(images, xy_position, ch, timepoint)
            
#             phase_labels = mint.get_otsu_mask(phase_image, min_area=20, max_area=5000)
            
#             for ch in fluor_images_dict:
#                 bkg_cor_image = bkg.back_sub(fluor_images_dict[ch], phase_labels>0, 
#                                              dilation=35, estimation_step=128, smoothing_sigma=60, show=False)
#                 fluor_images_dict[ch] = bkg_cor_image[0]
                
#             save_suffix = 'xy'+str(xy_position)+'_time'+str(timepoint)
#             mint.plot_cell_labels(phase_labels, save_suffix, save_path)
            
#             for cell_label in np.unique(phase_labels.ravel()):
#                 if cell_label > 0:
            
                
#                 oned_df, crop_pad, cropped_cell_mask = mint.get_medial_axis_for_label(phase_labels, lbl, radius_px=8, half_angle=25, 
#                                                                                  cap_knot=16, max_degree=30, verbose=True)
# # plot_intensity_profiles_and_cells(images, xy_position, timepoint, save_path)
            
    
    
    
