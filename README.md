# Martinez-Sarmiento_et_al_2023
Software used in the research manuscript "Super-Resolution Imaging Uncovers Key Temporal Changes in Chromatin Structure and Pluripotent Gene Reactivation in Single Cells Undergoing Reprogramming"

Files:

1- callVoronoiMultiSeg - MATLAB script to perform Voronoi Segmentation Analysis

2- PostVA_mod - MATLAB script that perform simple calculations on the VA analyzed data over multiple files.

3- Voronoi_batch analysis folder - Open the MATLAB file "to run" to run the code. It performs Voronoi Clustering-based outlier removal to filter STORM data automatically over multiple files. It performs Voronoi Tessellation to obtain the Voronoi areas and then perform Voronoi Clustering. From the clustered image, it removes outlier data that deviates from the normal distribution of nuclear area vs localization density. Finally, it combines the filtered data with the non-clustered data and run Voronoi segmentation again.

4- oniRawBincropper - MATLAB scrip that contain a Function to crop the a tiff file and overlay the corresponding bin file. It outputs a .Dax and a .bin file. Outputs a file for each cropped region

5- Stack_tif - MATLAB script that automatically stacks TIF files.

6- STORM Data Analysis Software -Lakadamyali lab- April 23, 2021 folder -It contains MATLAB scripts to perfom multiple calculation over single-molecule localization microscopy data. Open the files "data_analysis_software" to run the code.
