function [features_grid] = crossed_features_3D_grid(grid, all_crossed)
%CROSSED_FEATURES_3D_GRID computes the aggregated crossed values on the
%grid

% crossed features initialization
features_grid.centers = zeros(grid.nb_x,grid.nb_y,grid.nb_z,3);

features_grid.cells_number.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.cells_number.min = 0;
features_grid.cells_number.max = 0;
features_grid.cells_number.title = 'Cells Number';
features_grid.cells_number.scaling_factor = 1;

features_grid.corresp.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.corresp.min = 0;
features_grid.corresp.max = 0;
features_grid.corresp.title = 'Proportion Corresponding';
features_grid.corresp.scaling_factor = 1000;

features_grid.n_c_ratio.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.n_c_ratio.min = 0;
features_grid.n_c_ratio.max = 0;
features_grid.n_c_ratio.title = 'N/C Ratio';
features_grid.n_c_ratio.scaling_factor = 1000;

features_grid.centroid_shift.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.centroid_shift.min = 0;
features_grid.centroid_shift.max = 0;
features_grid.centroid_shift.title = 'Centroid Shift';
features_grid.centroid_shift.scaling_factor = 1;

features_grid.PC1_alignment.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.PC1_alignment.min = 0;
features_grid.PC1_alignment.max = 0;
features_grid.PC1_alignment.title = 'PC1 Alignment';
features_grid.PC1_alignment.scaling_factor = 1000;

features_grid.PC2_alignment.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.PC2_alignment.min = 0;
features_grid.PC2_alignment.max = 0;
features_grid.PC2_alignment.title = 'PC2 Alignment';
features_grid.PC2_alignment.scaling_factor = 1000;

features_grid.PC3_alignment.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.PC3_alignment.min = 0;
features_grid.PC3_alignment.max = 0;
features_grid.PC3_alignment.title = 'PC3 Alignment';
features_grid.PC3_alignment.scaling_factor = 1000;

% populating individual crossed features on the bins
for i = 1:grid.nb_x,
    for j = 1:grid.nb_y,
        for k = 1:grid.nb_z,
            indtemp1 = find(all_crossed(:,2)>=grid.x_bins(i) & all_crossed(:,2)<grid.x_bins(i+1) & ...
                all_crossed(:,3) >=grid.y_bins(j) & all_crossed(:,3)<grid.y_bins(j+1) & ...
                all_crossed(:,4)>=grid.z_bins(k) & all_crossed(:,4)<grid.z_bins(k+1));
            features_grid.cells_number.vals(i,j,k) = length(indtemp1);
            
            if features_grid.cells_number.vals(i,j,k)>= 10,
                features_grid.centers(i,j,k,:) = [mean(all_crossed(indtemp1,2)),mean(all_crossed(indtemp1,3)),mean(all_crossed(indtemp1,4))];
                features_grid.corresp.vals(i,j,k) = mean(all_crossed(indtemp1,5));
                features_grid.n_c_ratio.vals(i,j,k) = mean(all_crossed(indtemp1,6));
                features_grid.centroid_shift.vals(i,j,k) = mean(all_crossed(indtemp1,7));
                
                features_grid.PC1_alignment.vals(i,j,k) = mean(abs(dot(all_crossed(indtemp1,8:10),all_crossed(indtemp1,20:22),2)));
                features_grid.PC2_alignment.vals(i,j,k) = mean(abs(dot(all_crossed(indtemp1,11:13),all_crossed(indtemp1,23:25),2)));
                features_grid.PC3_alignment.vals(i,j,k) = mean(abs(dot(all_crossed(indtemp1,14:16),all_crossed(indtemp1,26:28),2)));
                
                
            end
        end
    end
end

vals_temp = features_grid.cells_number.vals(:);
indtemp = find(vals_temp>=10);


names = fieldnames(features_grid);

for i = 2:8,
    vals_temp1 = features_grid.(names{i}).vals;
    vals_temp = vals_temp1(:);
    features_grid.(names{i}).min = min(vals_temp(indtemp));
    features_grid.(names{i}).max = max(vals_temp(indtemp));
end    

end

