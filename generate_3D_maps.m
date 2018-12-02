function [] = generate_3D_maps(opt)
%GENERATE_3D_MAPS Generate 3D maps of nuclei and cells features
%   the input opt contains information about the type of maps that should
%   be generated.
%   the output of the function is a set of maps

gp = length(opt.path);
for i = 1:gp
    % compute features on the grid
    if opt.nuclei
        if exist([opt.path{i},'nuclei_grid.mat'], 'file') == 0
            % compute and save nuclei features on the 3D grid
            nuclei_features(opt,i);
        end
    end
    
    if opt.cells
        if exist([opt.path{i},'cells_grid.mat'], 'file') == 0
            % compute and save cell features on the 3D grid
            cells_features(opt,i);
        end
    end
    
    if opt.crossed
        if exist([opt.path{i},'crossed_grid.mat'], 'file') == 0
            % compute and save crossed features on the 3D grid
            crossed_features(opt,i);
        end
    end
end

% compute visualizations
% compute min and max
[visualization_nuclei,visualization_cells,visualization] = compute_min_max(opt);

% draw maps
for i = 1:gp
    if opt.nuclei
        % load nuclei features on the 3D grid
        res_saved = load([opt.path{i},'nuclei_grid.mat']);
        
        opt.save_folder = [opt.path{i},'figures/nuclei/3D_Maps/'];
        draw_maps_nuclei_cells(res_saved.nuclei_grid,res_saved.grid,visualization_nuclei,opt);
        
        opt.save_folder = [opt.path{i},'figures/nuclei/spatial_profiles/'];
        draw_spatial_profiles(res_saved.spatial_profiles,res_saved.grid_sp,opt);
    end
    
    if opt.cells
        % load cell features on the 3D grid
        res_saved = load([opt.path{i},'cells_grid.mat']);
        
        opt.save_folder = [opt.path{i},'figures/cells/3D_Maps/'];
        draw_maps_nuclei_cells(res_saved.cells_grid,res_saved.grid,visualization_cells,opt);
        
        opt.save_folder = [opt.path{i},'figures/cells/spatial_profiles/'];
        draw_spatial_profiles(res_saved.spatial_profiles,res_saved.grid_sp,opt);
    end
    
    if opt.crossed
        % load crossed features on the 3D grid
        res_saved = load([opt.path{i},'crossed_grid.mat']);
        
        opt.save_folder = [opt.path{i},'figures/crossed/3D_Maps/'];
        draw_maps_crossed(res_saved.crossed_grid,res_saved.grid,visualization,opt);
        
        opt.save_folder = [opt.path{i},'figures/crossed/spatial_profiles/'];
        draw_spatial_profiles(res_saved.spatial_profiles,res_saved.grid_sp,opt);
    end
end

end

