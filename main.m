% this script generates individuals and grid scale cells features and
% produces 3D maps and spatial profiles

%% add path


% add path
root_dir = '';

% Subfolders
addpath(fullfile(root_dir, 'data'));
addpath(fullfile(root_dir, 'utils'));

presets


%% several growth plates at a time

% write here the path to the folders containing the datasets for each
% growth plate
opt.path = {'data/Nuclei_and_Cells_DT_S18_m6_wt/';'data/Nuclei_and_Cells_PT_S18_m6_wt/'}

% this option enables to flip the x axis (1: flipped, 0: not flipped) -
% match the order of opt.path
opt.flip_x_axis = {0;1};

% dimensions of the grid
opt.delta_x = 75;
opt.delta_y = 75;
opt.delta_z = 75;

% 3D view
opt.az = -162;
opt.el = 59;
% maps for the nuclei - value: 0 or 1 (if they should be included)
opt.nuclei = 1;

% maps for the cells - value: 0 or 1 (if they should be included)
opt.cells = 1;

% maps for crossed features (alignment, shift etc..)
opt.crossed = 1;

% save figs
opt.save_figs = 1;

generate_3D_maps(opt);

%% DT and PT together

% write here the path to the folders containing the datasets for each
% growth plate
% NB: keep the corresponding order for DT and PT
opt.path_DT = {'data/Nuclei_and_Cells_DT_S18_m6_wt/';'data/Nuclei_and_Cells_DT_S18_m6_wt/';}
opt.path_PT = {'data/Nuclei_and_Cells_PT_S18_m6_wt/';'data/Nuclei_and_Cells_PT_S18_m6_wt/';}
opt.specimen_name = {'S18_m6_wt';'S18_m6_wt_bis';}

% this option enables to flip the x axis (1: flipped, 0: not flipped) -
% match the order of opt.path_DT 
opt.flip_x_axis_DT_PT = {0;1};

% dimensions of the grid
opt.delta_x = 75;
opt.delta_y = 75;
opt.delta_z = 75;

% 3D view
opt.az = -162;
opt.el = 59;

% adjustment DT for visualization purposes
opt.DT_shift_x = 1100;
opt.DT_shift_y = 150;

% maps for the nuclei - value: 0 or 1 (if they should be included)
opt.nuclei = 1;

% maps for the cells - value: 0 or 1 (if they should be included)
opt.cells = 1;

% maps for crossed features (alignment, shift etc..)
opt.crossed = 1;

% save figs
opt.save_figs = 1;

generate_3D_maps_DT_and_PT(opt);
