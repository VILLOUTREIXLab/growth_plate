function [visualization_nuclei,visualization_cells,visualization] = compute_min_max(opt)
%COMPUTE_MIN_MAX This function computes the range of values used to draw
%each maps
% find min and max

visualization.az = opt.az; visualization.el = opt.el;
visualization_nuclei = visualization;
visualization_cells = visualization;

gp = length(opt.path);

for k = 1:gp
    if opt.nuclei & opt.cells
        
        % min max is computed accross the two modalities
        % load nuclei features on the 3D grid
        res_saved = load([opt.path{k},'nuclei_grid.mat']);
        nuclei_grid = res_saved.nuclei_grid;
        grid = res_saved.grid;
        clear res_saved;
        
        names = fieldnames(nuclei_grid);
        for i = 2:length(names),
            if k == 1,
                visualization_nuclei.(names{i}).min = nuclei_grid.(names{i}).min;
                visualization_nuclei.(names{i}).max = nuclei_grid.(names{i}).max;
            else
                visualization_nuclei.(names{i}).min = min(visualization_nuclei.(names{i}).min,nuclei_grid.(names{i}).min);
                visualization_nuclei.(names{i}).max = max(visualization_nuclei.(names{i}).max,nuclei_grid.(names{i}).max);
            end
        end
        
        % load cell features on the 3D grid
        res_saved = load([opt.path{k},'cells_grid.mat']);
        cells_grid = res_saved.cells_grid;
        grid = res_saved.grid;
        clear res_saved;
        
        names = fieldnames(cells_grid);
        for i = 2:length(names),
            if k == 1,
                visualization_cells.(names{i}).min = cells_grid.(names{i}).min;
                visualization_cells.(names{i}).max = cells_grid.(names{i}).max;
            else
                visualization_cells.(names{i}).min = min(visualization_cells.(names{i}).min,cells_grid.(names{i}).min);
                visualization_cells.(names{i}).max = max(visualization_cells.(names{i}).max,cells_grid.(names{i}).max);
            end
        end
        i = 5;
        visualization_cells.(names{i}).min = min(visualization_nuclei.(names{i}).min,cells_grid.(names{i}).min);
        visualization_nuclei.(names{i}).min = visualization_cells.(names{i}).min;
        visualization_cells.(names{i}).max = max(visualization_nuclei.(names{i}).max,cells_grid.(names{i}).max);
        visualization_nuclei.(names{i}).max = visualization_cells.(names{i}).max;
        
        for i = length(names)-3:length(names),
            visualization_cells.(names{i}).min = min(visualization_nuclei.(names{i}).min,cells_grid.(names{i}).min);
            visualization_nuclei.(names{i}).min = visualization_cells.(names{i}).min;
            visualization_cells.(names{i}).max = max(visualization_nuclei.(names{i}).max,cells_grid.(names{i}).max);
            visualization_nuclei.(names{i}).max = visualization_cells.(names{i}).max;
        end
        
    elseif opt.nuclei & (opt.cells == 0)
        
        % load nuclei features on the 3D grid
        res_saved = load([opt.path{k},'nuclei_grid.mat']);
        nuclei_grid = res_saved.nuclei_grid;
        grid = res_saved.grid;
        clear res_saved;
        
        names = fieldnames(nuclei_grid);
        for i = 2:length(names),
            if k == 1
                visualization_nuclei.(names{i}).min = nuclei_grid.(names{i}).min;
                visualization_nuclei.(names{i}).max = nuclei_grid.(names{i}).max;
            else
                visualization_nuclei.(names{i}).min = min(visualization_nuclei.(names{i}).min,nuclei_grid.(names{i}).min);
                visualization_nuclei.(names{i}).max = max(visualization_nuclei.(names{i}).max,nuclei_grid.(names{i}).max);
            end
        end
        
    elseif (opt.nuclei == 0) & opt.cells
        
        % load cell features on the 3D grid
        res_saved = load([opt.path{k},'cells_grid.mat']);
        cells_grid = res_saved.cells_grid;
        grid = res_saved.grid;
        clear res_saved;
        
        names = fieldnames(cells_grid);
        for i = 2:length(names),
            if k == 1,
                visualization_cells.(names{i}).min = cells_grid.(names{i}).min;
                visualization_cells.(names{i}).max = cells_grid.(names{i}).max;
            else
                visualization_cells.(names{i}).min = min(visualization_cells.(names{i}).min,cells_grid.(names{i}).min);
                visualization_cells.(names{i}).max = max(visualization_cells.(names{i}).max,cells_grid.(names{i}).max);
            end
        end
    end
    
    
    if opt.crossed
        % load crossed features on the 3D grid
        res_saved = load([opt.path{k},'crossed_grid.mat']);
        crossed_grid = res_saved.crossed_grid;
        grid = res_saved.grid;
        clear res_saved;
        
        names = fieldnames(crossed_grid);
        for i = 2:(length(names)-3),
            if k == 1
                visualization.(names{i}).min = crossed_grid.(names{i}).min;
                visualization.(names{i}).max = crossed_grid.(names{i}).max;
            else
                visualization.(names{i}).min = min(visualization.(names{i}).min,crossed_grid.(names{i}).min);
                visualization.(names{i}).max = max(visualization.(names{i}).max,crossed_grid.(names{i}).max);
            end
        end
        
        temp_min = [crossed_grid.(names{(length(names)-2)}).min,crossed_grid.(names{(length(names)-1)}).min,crossed_grid.(names{length(names)}).min];
        temp_max = [crossed_grid.(names{(length(names)-2)}).max,crossed_grid.(names{(length(names)-1)}).max,crossed_grid.(names{length(names)}).max];
        
        for i = (length(names)-2):length(names),
            if k == 1
                visualization.(names{i}).min = min(temp_min);
                visualization.(names{i}).max = max(temp_max);
            else
                visualization.(names{i}).min = min(visualization.(names{i}).min,min(temp_min));
                visualization.(names{i}).max = max(visualization.(names{i}).max,max(temp_max));
            end
        end
    end
end

end

