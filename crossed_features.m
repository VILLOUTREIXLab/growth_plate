function [] = crossed_features(opt,gp)
%CROSSED_FEATURES Computes individual crossed features and their aggregated
%values on a 3D grid

if exist([opt.path{gp},'all_crossed.mat'], 'file') == 0
    
    % load coordinates
    [numbers,txt,raw] = xlsread([opt.path{gp},'Tile_coordinates.xlsx']);
    
    coordinates = zeros(size(txt,1)-3,5);
    for i = 4:size(txt,1),
        temp =  char(txt(i,1));
        res = strsplit(temp,'_POS');
        coordinates(i-3,1) = str2num(char(res(2)));
        coordinates(i-3,2:5) = numbers(i-3,:);
    end
    
    disp('Computing individual crossed features')
    
    % compute individual cell features
    all_crossed = [];
    
    % loop over all the positions
    for position = coordinates(:,1)'
        disp(num2str(position));
        
        % load individual cell features
        load([opt.path{gp},'c_n_pos',num2str(position),' (Characteristics).mat']);
        
        % here need to check if nuclear cell ratio is good
        indtemp = find(G.inter.volume_ratio>1);
        if length(indtemp)>0,
            
            all_crossed_temp = zeros(length(indtemp),7);
            all_crossed_temp(:,1) = position*ones(length(indtemp),1);
            all_crossed_temp(:,2:4) = G.cel.centroids(indtemp,:);
            
            % organizing individual cell position in the global coordinate
            % system
            indtemp1 = find(coordinates(:,1) == position);
            all_crossed_temp(:,2) = all_crossed_temp(:,2)+coordinates(indtemp1,3)*ones(size(all_crossed_temp,1),1);
            all_crossed_temp(:,3) = all_crossed_temp(:,3)-coordinates(indtemp1,2)*ones(size(all_crossed_temp,1),1);
            all_crossed_temp(:,4) = all_crossed_temp(:,4)+coordinates(indtemp1,4)*ones(size(all_crossed_temp,1),1);
            
            % proportion of corresponding cells compared to number of
            % nuclei only
            indtemp1 = find(G.nuc.index);
            all_crossed_temp(:,5) = (length(indtemp)/length(indtemp1))*ones(length(indtemp),1);
            
            % n/c ratio
            all_crossed_temp(:,6) = 1./G.inter.volume_ratio(indtemp);
                        
            % Centroid shift
            all_crossed_temp(:,7) = sqrt(G.inter.centroid_shift_physical_coords(indtemp,1).^2+G.inter.centroid_shift_physical_coords(indtemp,2).^2+G.inter.centroid_shift_physical_coords(indtemp,3).^2);
            
            % PC nuclei
            vals = cell(1,3);
            for pc = 1 : 3
                tmp = cell2mat(cellfun(@(x) x(:,pc)', G.nuc.PCA_coeff(indtemp), 'UniformOutput', false));
                vals{pc} = [vals{pc}; tmp];
            end
            
            all_crossed_temp(:,8:10) = vals{1}(:,1:3);
            all_crossed_temp(:,11:13) = vals{2}(:,1:3);
            all_crossed_temp(:,14:16) = vals{3}(:,1:3);
            all_crossed_temp(:,17:19) = G.nuc.PCA_latent(indtemp,:);
            
            
            % PC Cells
            vals = cell(1,3);
            for pc = 1 : 3
                tmp = cell2mat(cellfun(@(x) x(:,pc)', G.cel.PCA_coeff(indtemp), 'UniformOutput', false));
                vals{pc} = [vals{pc}; tmp];
            end
            
            all_crossed_temp(:,20:22) = vals{1}(:,1:3);
            all_crossed_temp(:,23:25) = vals{2}(:,1:3);
            all_crossed_temp(:,26:28) = vals{3}(:,1:3);
            all_crossed_temp(:,29:31) = G.cel.PCA_latent(indtemp,:);
            
            all_crossed = [all_crossed;all_crossed_temp];
        end
    end
    
    disp('Done');
    
    save([opt.path{gp},'all_crossed.mat'],'all_crossed');
else
    load([opt.path{gp},'all_crossed.mat']);
end

% compute crossed features on 3D grid
characteristics_on_grid_crossed(all_crossed,opt,gp);

end