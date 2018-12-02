function [features_grid] = features_3D_grid(grid,all,type)
%FEATURES_3D_GRID aggregates the individual cell features on the 3D grid
%   Detailed explanation goes here

% nuclei features initialization
features_grid.centers = zeros(grid.nb_x,grid.nb_y,grid.nb_z,3);

features_grid.object_number.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.object_number.min = 0;
features_grid.object_number.max = 0;
if type == 1,
    features_grid.object_number.title = 'Nuclei Number';
elseif type == 2,
    features_grid.object_number.title = 'Cell Number';
end
features_grid.object_number.scaling_factor = 1;

features_grid.volume.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.volume.min = 0;
features_grid.volume.max = 0;
features_grid.volume.title = 'Volume';
features_grid.volume.scaling_factor = 1;
features_grid.volume.scaling_factor_log = 100;

features_grid.surface_area.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.surface_area.min = 0;
features_grid.surface_area.max = 0;
features_grid.surface_area.title = 'Surface Area';
features_grid.surface_area.scaling_factor = 1;
features_grid.surface_area.scaling_factor_log = 100;

features_grid.sphericity.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.sphericity.min = 0;
features_grid.sphericity.max = 0;
features_grid.sphericity.title = 'Sphericity';
features_grid.sphericity.scaling_factor = 10000;

features_grid.occupation.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.occupation.min = 0;
features_grid.occupation.max = 0;
features_grid.occupation.title = 'Occupation';
features_grid.occupation.scaling_factor = 10000;

features_grid.density.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z);
features_grid.density.min = 0;
features_grid.density.max = 0;
features_grid.density.title = 'Density';
features_grid.density.scaling_factor = 10^8;

features_grid.density_delaunay.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.density_delaunay.min = 0;
features_grid.density_delaunay.max = 0;
features_grid.density_delaunay.title = 'Delaunay Density';
features_grid.density_delaunay.scaling_factor = 10^8;

features_grid.PC1.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.PC1.min = 0;
features_grid.PC1.max = 0;
features_grid.PC1.title = 'PC1 coefficient';
features_grid.PC1.scaling_factor = 10^2;
features_grid.PC1.orientations = zeros(grid.nb_x,grid.nb_y,grid.nb_z,3);
features_grid.PC1.sph_var = zeros(grid.nb_x,grid.nb_y,grid.nb_z);

features_grid.PC2.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.PC2.min = 0;
features_grid.PC2.max = 0;
features_grid.PC2.title = 'PC2 coefficient';
features_grid.PC2.scaling_factor = 10^2;
features_grid.PC2.orientations = zeros(grid.nb_x,grid.nb_y,grid.nb_z,3);
features_grid.PC2.sph_var = zeros(grid.nb_x,grid.nb_y,grid.nb_z);

features_grid.PC3.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.PC3.min = 0;
features_grid.PC3.max = 0;
features_grid.PC3.title = 'PC3 coefficient';
features_grid.PC3.scaling_factor = 10^2;
features_grid.PC3.orientations = zeros(grid.nb_x,grid.nb_y,grid.nb_z,3);
features_grid.PC3.sph_var = zeros(grid.nb_x,grid.nb_y,grid.nb_z);

features_grid.ratio_PC1_PC2.vals = zeros(grid.nb_x,grid.nb_y,grid.nb_z,2);
features_grid.ratio_PC1_PC2.min = 0;
features_grid.ratio_PC1_PC2.max = 0;
features_grid.ratio_PC1_PC2.title = 'Ratio PC2/PC1';
features_grid.ratio_PC1_PC2.scaling_factor = 10^3;

% populating individual cell features on the bin
for i = 1:grid.nb_x,
    for j = 1:grid.nb_y,
        for k = 1:grid.nb_z,
            indtemp1 = find( all(:,5)>=grid.x_bins(i) & all(:,5)<grid.x_bins(i+1) & ...
                all(:,6) >=grid.y_bins(j) & all(:,6)<grid.y_bins(j+1) & ...
                all(:,7)>=grid.z_bins(k) & all(:,7)<grid.z_bins(k+1));
            features_grid.object_number.vals(i,j,k) = length(indtemp1);
            
            if features_grid.object_number.vals(i,j,k)>= 10,
                features_grid.centers(i,j,k,:) = [mean(all(indtemp1,5)),mean(all(indtemp1,6)),mean(all(indtemp1,7))];
                features_grid.volume.vals(i,j,k,1) = mean(all(indtemp1,2));
                features_grid.volume.vals(i,j,k,2) = std(all(indtemp1,2));
                features_grid.surface_area.vals(i,j,k,1) = mean(all(indtemp1,3));
                features_grid.surface_area.vals(i,j,k,2) = std(all(indtemp1,3));
                features_grid.sphericity.vals(i,j,k,1) = mean(all(indtemp1,4));
                features_grid.sphericity.vals(i,j,k,2) = std(all(indtemp1,4));
                features_grid.occupation.vals(i,j,k,1) = sum(all(indtemp1,2))/grid.volume_bin;
                features_grid.density.vals(i,j,k) = features_grid.object_number.vals(i,j,k)/grid.volume_bin;
                
                features_grid.density_delaunay.vals(i,j,k,1) = mean(all(indtemp1,20));
                features_grid.density_delaunay.vals(i,j,k,2) = std(all(indtemp1,20));
                
                features_grid.PC1.vals(i,j,k,1) = mean(all(indtemp1,17));
                features_grid.PC1.vals(i,j,k,2) = std(all(indtemp1,17));
                features_grid.PC2.vals(i,j,k,1) = mean(all(indtemp1,18));
                features_grid.PC2.vals(i,j,k,2) = std(all(indtemp1,18));
                features_grid.PC3.vals(i,j,k,1) = mean(all(indtemp1,19));
                features_grid.PC3.vals(i,j,k,2) = std(all(indtemp1,19));
                
                features_grid.ratio_PC1_PC2.vals(i,j,k,1) = mean(all(indtemp1,18)./all(indtemp1,17));
                features_grid.ratio_PC1_PC2.vals(i,j,k,1) = std(all(indtemp1,18)./all(indtemp1,17));
                
                % Orientations
                % PC1
                vals = all(indtemp1,8:10);
                
                % making sure that each orientation is represented in both directions:
                vals = [vals; -vals];
                
                % approximating the average orientation of the nuclei as the PC1 of all orientations (i.e. the PC1 of all the PC1's of the nuclei ? and their negative):
                [coeff,score,latent] = pca(vals);
                
                % finding the orientations in the same hemisphere as the PC1 of all the
                % orientations
                indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
                
                % spherical mean
                R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
                
                features_grid.PC1.orientations(i,j,k,:) = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];
                
                % spherical variance
                features_grid.PC1.sph_var(i,j,k) = (length(indtemp2) - R)/length(indtemp2);
                
                % PC2
                vals = all(indtemp1,11:13);
                
                % making sure that each orientation is represented in both directions:
                vals = [vals; -vals];
                
                % approximating the average orientation of the nuclei as the PC1 of all orientations (i.e. the PC1 of all the PC2's of the nuclei and their negative):
                [coeff,score,latent] = pca(vals);
                
                % finding the orientations in the same hemisphere as the PC1 of all the
                % orientations
                indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
                
                % spherical mean
                R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
                
                features_grid.PC2.orientations(i,j,k,:) = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];
                
                % spherical variance
                features_grid.PC2.sph_var(i,j,k) = (length(indtemp2) - R)/length(indtemp2);
                
                
                % PC3
                vals = all(indtemp1,14:16);
                
                % making sure that each orientation is represented in both directions:
                vals = [vals; -vals];
                
                % approximating the average orientation of the nuclei as the PC1 of all orientations (i.e. the PC1 of all the PC3's of the nuclei and their negative):
                [coeff,score,latent] = pca(vals);
                
                % finding the orientations in the same hemisphere as the PC1 of all the
                % orientations
                indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
                
                % spherical mean
                R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
                features_grid.PC3.orientations(i,j,k,:) = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];
                
                % spherical variance
                features_grid.PC3.sph_var(i,j,k) = (length(indtemp2) - R)/length(indtemp2);
                
            end
        end
    end
end

vals_temp = features_grid.object_number.vals(:);
indtemp = find(vals_temp>=10);

names = fieldnames(features_grid);

for i = 2:length(names),
    vals_temp1 = features_grid.(names{i}).vals;
    vals_temp = vals_temp1(:);
    features_grid.(names{i}).min = min(vals_temp(indtemp))
    features_grid.(names{i}).max = max(vals_temp(indtemp))
end    
end

