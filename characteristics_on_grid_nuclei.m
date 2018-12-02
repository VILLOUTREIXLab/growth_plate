function [] = characteristics_on_grid_nuclei(all_cells_nuclei,opt,gp)
%CHARACTERISTICS_ON_GRID_NUCLEI Computes aggregated nuclei features on a 3D
%grid
%   Input: all individual nuclei features
%   Output: aggregated nuclei features on 3D grid

% each row of all_cells_nuclei is nucleus
% the columns contains the individual nuclei features
% 1 - stack id
% 2 - volume
% 3 - surface area
% 4 - sphericity
% 5-7 - centroid x,y,z coordinates
% 8-10 - PC1 x,y,z orientation
% 11-13 - PC2 x,y,z orientation
% 14-16 - PC3 x,y,z orientation
% 17-19 - PC1,PC2,PC3 latent coefficient
% 20 - Delaunay density

% computing grid features
grid.min_x = min(all_cells_nuclei(:,5));
grid.max_x = max(all_cells_nuclei(:,5));
grid.nb_x = round((grid.max_x-grid.min_x)/opt.delta_x);
grid.x_bins = grid.min_x:opt.delta_x:(grid.min_x+grid.nb_x*opt.delta_x);
grid.min_y = min(all_cells_nuclei(:,6));
grid.max_y = max(all_cells_nuclei(:,6));
grid.nb_y = round((grid.max_y-grid.min_y)/opt.delta_y);
grid.y_bins = grid.min_y:opt.delta_y:(grid.min_y+grid.nb_y*opt.delta_y);
grid.min_z = min(all_cells_nuclei(:,7));
grid.max_z = max(all_cells_nuclei(:,7));
grid.nb_z = round((grid.max_z-grid.min_z)/opt.delta_z);
grid.z_bins = grid.min_z:opt.delta_z:(grid.min_z+grid.nb_z*opt.delta_z);
grid.volume_bin = opt.delta_x*opt.delta_y*opt.delta_z;

nuclei_grid = features_3D_grid(grid,all_cells_nuclei,1);

% computing grid features for spatial profile
grid_sp.min_x = min(all_cells_nuclei(:,5));
grid_sp.max_x = max(all_cells_nuclei(:,5));
grid_sp.nb_x = round((grid_sp.max_x-grid_sp.min_x)/opt.delta_x);
grid_sp.x_bins = grid_sp.min_x:opt.delta_x:(grid_sp.min_x+grid_sp.nb_x*opt.delta_x);
% define 3 zones for y
grid_sp.min_y = min(all_cells_nuclei(:,6));
grid_sp.max_y = max(all_cells_nuclei(:,6));
grid_sp.nb_y = 3;
grid_sp.y_bins = grid_sp.min_y:(grid_sp.max_y-grid_sp.min_y)/3:grid_sp.max_y;
% define 3 zones for z
grid_sp.min_z = min(all_cells_nuclei(:,7));
grid_sp.max_z = max(all_cells_nuclei(:,7));
grid_sp.nb_z = 3;
grid_sp.z_bins = grid_sp.min_z:(grid_sp.max_z-grid_sp.min_z)/3:grid_sp.max_z;
grid_sp.volume_bin = opt.delta_x*opt.delta_y*opt.delta_z;

spatial_profiles = features_3D_grid(grid_sp,all_cells_nuclei,2);

save([opt.path{gp},'nuclei_grid.mat'],'nuclei_grid','grid','spatial_profiles','grid_sp');

end

