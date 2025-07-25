# PyntensityProfile

Author: Alexandros Papagiannakis, HHMI @Stanford University, Christine Jacobs-Wagner lab, Sarafan ChEM-H, 2025

 Cite:
    
    Nonequilibrium polysome dynamics promote chromosome segregation and its coupling to cell growth in Escherichia coli
    
    Alexandros Papagiannakis, Qiwei Yu, Sander K. Govers, Wei-Hsiang Lin,  Ned S. Wingreen, Christine Jacobs-Wagner
    
    eLife2025;14:RP104276 DOI: https://doi.org/10.7554/eLife.104276.3

## Instructions
This repository includes functions that can be used to generate intensity profiles from single cells.

In otder to use this package first load the nd2 images using the function load_images() for a specified image directory.
This package works for nd2 images with the 'mct' iteration axis, which correspond to a timelapse movie with multiple xy positions and channels.

The loaded images can be stored to a viariable called images. Then for a specified xy_position and timepoint the function below should be executed:
plot_intensity_profiles_and_cells(images, xy_position, timepoint, save_path)

The save_path parameter is the path to the directory where the cell images and intensity profile will be stored.

<br> Once the plot_intensity_profiles_and_cells() function is run the user will have the opportubity to select the cell label for which they want to generate the intensity profile: </br>
<img src="https://github.com/user-attachments/assets/623eb748-ade5-4b8f-bc2d-668bf1bcdee1" align="center" width="500"/>

<br>By selecting an integer label as input:</br>
<img src="https://github.com/user-attachments/assets/21c29d9a-e27d-436b-87a4-a5a3d6592b16" align="center" width="500"/>

<br>The medial axis will be drawn for the selected label:</br>
<img src="https://github.com/user-attachments/assets/4fcfb3a9-348d-40f2-9f59-d9ec19212626" align="center" width="100"/>

<br>And the intensity profiles will be generated for the fluorescent channels:</br>
<img src="https://github.com/user-attachments/assets/9283aab3-4899-4551-b878-6d52a5c474f7" align="center" width="400"/>




