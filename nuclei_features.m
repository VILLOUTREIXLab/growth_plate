function [] = nuclei_features(opt,gp)
%NUCLEI_FEATURES Computes individual nuclei features and their aggregated
%values on a 3D grid

if exist([opt.path{gp},'all_cells_nuclei.mat'], 'file') == 0
    
    % load coordinates
    [numbers,txt,raw] = xlsread([opt.path{gp},'Tile_coordinates.xlsx']);
    
    coordinates = zeros(size(txt,1)-3,5);
    for i = 4:size(txt,1),
        temp =  char(txt(i,1));
        res = strsplit(temp,'_POS');
        coordinates(i-3,1) = str2num(char(res(2)));
        coordinates(i-3,2:5) = numbers(i-3,:);
    end
    
    disp('Computing individual nuclei features')
    
    % compute individual nuclei features
    all_cells_nuclei = [];
    
    % loop over all the positions
    for position = coordinates(:,1)'
        disp(num2str(position));
        
        % load individual nuclei features
        load([opt.path{gp},'c_n_pos',num2str(position),' (Characteristics).mat']);
        
        indtemp = find(G.nuc.index);
        if length(indtemp)>0,
            
            all_cells_nuclei_temp = zeros(length(indtemp),7);
            all_cells_nuclei_temp(:,1) = position*ones(length(indtemp),1);
            all_cells_nuclei_temp(:,2:4) = [G.nuc.volume(indtemp),G.nuc.surface_area(indtemp),G.nuc.sphericity(indtemp)];
            all_cells_nuclei_temp(:,5:7) = G.nuc.centroids(indtemp,:);
            
            vals = cell(1,3);
            for pc = 1 : 3
                tmp = cell2mat(cellfun(@(x) x(:,pc)', G.nuc.PCA_coeff(indtemp), 'UniformOutput', false));
                vals{pc} = [vals{pc}; tmp];
            end
            
            all_cells_nuclei_temp(:,8:10) = vals{1}(:,1:3);
            all_cells_nuclei_temp(:,11:13) = vals{2}(:,1:3);
            all_cells_nuclei_temp(:,14:16) = vals{3}(:,1:3);
            all_cells_nuclei_temp(:,17:19) = G.nuc.PCA_latent(indtemp,:);
            
            % organizing individual cell position in the global coordinate
            % system
            indtemp1 = find(coordinates(:,1) == position);
            
            all_cells_nuclei_temp(:,5) = -all_cells_nuclei_temp(:,5)+coordinates(indtemp1,3)*ones(size(all_cells_nuclei_temp,1),1);
            all_cells_nuclei_temp(:,8)  = -all_cells_nuclei_temp(:,8);
            all_cells_nuclei_temp(:,11) = -all_cells_nuclei_temp(:,11);
            all_cells_nuclei_temp(:,14) = -all_cells_nuclei_temp(:,14);
            all_cells_nuclei_temp(:,6) = all_cells_nuclei_temp(:,6)-coordinates(indtemp1,2)*ones(size(all_cells_nuclei_temp,1),1);
            all_cells_nuclei_temp(:,7) = all_cells_nuclei_temp(:,7)+coordinates(indtemp1,4)*ones(size(all_cells_nuclei_temp,1),1);
            
            all_cells_nuclei = [all_cells_nuclei;all_cells_nuclei_temp];
        end
    end
    
    
    if exist([opt.path{gp},'point_density_nuclei.mat'], 'file') == 0
        disp('Computing Delaunay density ...')
        [V,X,Y,Z, point_density] = calculate_nuclei_density(all_cells_nuclei(:,[5,6,7]), [1, 1, 1], 2);
        save([opt.path{gp},'point_density_nuclei.mat'],'V','X','Y','Z','point_density');
    end
    
    res_pd = load([opt.path{gp},'point_density_nuclei.mat']);
    all_cells_nuclei(:,20) = res_pd.point_density;
    
    disp('Done');
    
    save([opt.path{gp},'all_cells_nuclei.mat'],'all_cells_nuclei');
else
    load([opt.path{gp},'all_cells_nuclei.mat']);
end

% compute nuclei features on 3D grid
characteristics_on_grid_nuclei(all_cells_nuclei,opt,gp);

end
