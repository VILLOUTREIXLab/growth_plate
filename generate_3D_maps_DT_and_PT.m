function [] = generate_3D_maps_DT_and_PT(opt)
%GENERATE_3D_MAPS_DT_and_PT Generate 3D maps of nuclei and cells features
%   the input opt contains information about the type of maps that should
%   be generated.
%   the output of the function is a set of maps

opt.path = [opt.path_DT;opt.path_PT];
gp = length(opt.path);
for i = 1:gp
    % compute features on the grid
    if opt.nuclei
        if exist([opt.path{i},'nuclei_grid.mat'], 'file') == 0
            % compute and save nuclei features on the 3D grid
            nuclei_features(opt,gp);
        end
    end
    
    if opt.cells
        if exist([opt.path{i},'cells_grid.mat'], 'file') == 0
            % compute and save cell features on the 3D grid
            cells_features(opt,gp);
        end
    end
    
    if opt.crossed
        if exist([opt.path{i},'crossed_grid.mat'], 'file') == 0
            % compute and save crossed features on the 3D grid
            crossed_features(opt,gp);
        end
    end
end

% compute visualizations
% compute min and max
[visualization_nuclei,visualization_cells,visualization] = compute_min_max(opt);

b = length(opt.path_DT);

% draw maps
for i = 1:b
    
    if opt.nuclei
        opt.save_folder = ['data/',opt.specimen_name{i},'/DT_and_PT/nuclei/3D_Maps/'];
        draw_maps_nuclei_cells_DT_and_PT(visualization_nuclei,opt,1,i);
    end
    
    if opt.cells
        opt.save_folder = ['data/',opt.specimen_name{i},'/DT_and_PT/cells/3D_Maps/'];
        draw_maps_nuclei_cells_DT_and_PT(visualization_cells,opt,2,i);
        
    end
    
    if opt.crossed
        opt.save_folder = ['data/',opt.specimen_name{i},'/DT_and_PT/crossed/3D_Maps/'];
        draw_maps_crossed_DT_and_PT(visualization,opt,i);
    end
end

end

