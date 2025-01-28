# Intensity_profile
This repository includes functions that can be used to generate intensity profiles from single cells.

In otder to use this package first load the nd2 images using the function load_images() for a specified image directory.
This package works for nd2 images with the 'mct' iteration axis, which correspond to a timelapse movie with multiple xy positions and channels.

The loaded images can be stored to a viariable called images. Then for a specified xy_position and timepoint the function below should be executed:
plot_intensity_profiles_and_cells(images, xy_position, timepoint, save_path)

The save_path parameter is the path to the directory where the cell images and intensity profile will be stored.

Once the plot_intensity_profiles_and_cells() function is run the user will have the opportubity to select the cell label for which they want to generate the intensity profile:

![image](https://github.com/user-attachments/assets/623eb748-ade5-4b8f-bc2d-668bf1bcdee1)

