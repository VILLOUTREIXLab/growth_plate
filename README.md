# growth_plate
This repository contains the Matlab scripts used to generate 3D Maps of the growth plates.


# How to run the scripts to generate the 3D maps

1 - be sure to store each of the datasets in a subfolder of the folder data, e.g. ```Nuclei_and_Cells_DT_S18_m6_wt```
  * the files should be named ```c_n_pos3 (Characteristics).mat```, for each of the positions
  * the subfolder should contain an excel file named ```Tile_coordinates.xlsx``` which provides the coordinates of each of the positions
  * the datasets are prepared (preprocessed) using the scripts contained in the subfolder ```prepare_data/```. The first script to use is ```prepare_segmentation_file_for_analysis.``` who's input file is a binary tif image (segmented nuclei or cells). The second script to use is ```calc_morphological_characteristics_03_n_c_ratios.m``` who's input is the ' (SURFACES).mat file produced from the previous script
  
2 - edit the file ```main.m``` with the following options
  * names of the paths to the subfolders, e.g. ``` opt.path = {'data/Nuclei_and_Cells_DT_S18_m6_wt/';'data/Nuclei_and_Cells_PT_S18_m6_wt/'}; ``` for the two growth plates "DT_S18_m6" and "PT_S18_m6"
  * visualization options (az, el, visualizing nuclei, cells, crossed features)
  * NB: when drawing combined maps for DT and PT, the names of the paths are stored in the following way ```opt.path_DT = {'data/Nuclei_and_Cells_DT_S18_m6_wt/';'data/Nuclei_and_Cells_DT_S18_m6_wt/';}``` and ```opt.path_PT = {'data/Nuclei_and_Cells_PT_S18_m6_wt/';'data/Nuclei_and_Cells_PT_S18_m6_wt/';}``` - the order of the names in the lists should be the same for PT and DT
  
3 - run the file ```main.m```

4 - For any question, please contact me: paul.villoutreix@gmail.com
