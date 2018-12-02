function [crossed_grid, grid] = characteristics_on_grid_crossed(all_crossed,opt,gp)
%CHARACTERISTICS_ON_GRID_CROSSED Computes aggregated crossed features 
%(between cells and nuclei) on a 3D grid
%   Input: all individual cells crossed features
%   Output: aggregated crossed2 features on 3D grid

% each row of all_crossed is nucleus
% the columns contains the individual nuclei features
% 1 - stack id
% 2-4 - centroid x,y,z coordinates
% 5 - proportion corresponding
% 6 - n/c ratio
% 7 - centroid shift
% 8-10 - PC1 x,y,z orientation - nucleus
% 11-13 - PC1 x,y,z orientation - nucleus
% 14-16 - PC1 x,y,z orientation - nucleus
% 17-19 - PC1,PC2,PC3 latent coefficient - nucleus
% 20-22 - PC1 x,y,z orientation - cells
% 23-25 - PC1 x,y,z orientation - cells
% 26-28 - PC1 x,y,z orientation - cells
% 29-31 - PC1,PC2,PC3 latent coefficient - cells

% computing grid features
grid.min_x = min(all_crossed(:,2));
grid.max_x = max(all_crossed(:,2));
grid.nb_x = round((grid.max_x-grid.min_x)/opt.delta_x);
grid.x_bins = grid.min_x:opt.delta_x:(grid.min_x+grid.nb_x*opt.delta_x);
grid.min_y = min(all_crossed(:,3));
grid.max_y = max(all_crossed(:,3));
grid.nb_y = round((grid.max_y-grid.min_y)/opt.delta_y);
grid.y_bins = grid.min_y:opt.delta_y:(grid.min_y+grid.nb_y*opt.delta_y);
grid.min_z = min(all_crossed(:,4));
grid.max_z = max(all_crossed(:,4));
grid.nb_z = round((grid.max_z-grid.min_z)/opt.delta_z);
grid.z_bins = grid.min_z:opt.delta_z:(grid.min_z+grid.nb_z*opt.delta_z);
grid.volume_bin = opt.delta_x*opt.delta_y*opt.delta_z;

crossed_grid = crossed_features_3D_grid(grid, all_crossed);

% computing grid features for spatial profile
grid_sp.min_x = min(all_crossed(:,2));
grid_sp.max_x = max(all_crossed(:,2));
grid_sp.nb_x = round((grid_sp.max_x-grid_sp.min_x)/opt.delta_x);
grid_sp.x_bins = grid_sp.min_x:opt.delta_x:(grid_sp.min_x+grid_sp.nb_x*opt.delta_x);
% define 3 zones for y
grid_sp.min_y = min(all_crossed(:,3));
grid_sp.max_y = max(all_crossed(:,3));
grid_sp.nb_y = 3;
grid_sp.y_bins = grid_sp.min_y:(grid_sp.max_y-grid_sp.min_y)/3:grid_sp.max_y;
% define 3 zones for z
grid_sp.min_z = min(all_crossed(:,4));
grid_sp.max_z = max(all_crossed(:,4));
grid_sp.nb_z = 3;
grid_sp.z_bins = grid_sp.min_z:(grid_sp.max_z-grid_sp.min_z)/3:grid_sp.max_z;
grid_sp.volume_bin = opt.delta_x*opt.delta_y*opt.delta_z;

spatial_profiles = crossed_features_3D_grid(grid_sp,all_crossed);

                
save([opt.path{gp},'crossed_grid.mat'],'crossed_grid','grid','spatial_profiles','grid_sp');
end

