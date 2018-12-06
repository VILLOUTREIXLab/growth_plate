function prepare_segmentation_file_for_analysis

%{
INSTRUCTIONS:
This function accepts as input the .tif file with the segmented nuclei / cells, isolates only nuclei / cells that are valid in terms of their size, not
overlapping with the borders of the image, etc. Then, it prepares the data for the next script: calc_morphological_characteristics.
The outputs of this function are:
1. A .vff file with a name that ends with (Screened), and contains all the cells/nuclei that were left after all the screening.
2. A .mat file and a .vtk file that end with (Surfaces). The .vtk one is the input of the next script.
You should run this script 4 times if you use only nuclei (once for each zone), or 8 times if you're also using cells.

Spacing info: hz = [0.091 0.091 0.654], phz, pz, rz = [0.091 0.091 0.470]
E16.5 samples = [0.194 0.194 0.387]
%}

sd_to_size_coeff = 2 * 2.354; % constant - do not change!
last_spacing_filename = 'last_spacing.mat'; % constant - do not change!
n_parallel_workers = 22; % to get info about the availability of workers type 'get_num_of_workers', in the command window.
opengl_mode = 'Hardware';
min_volume = 70; % in cubic microns rz = 300 pz = 300  hz = 1000
max_volume = 5000; % rz = 800 PZ = 1000 HZ = 9000 in cubic microns % if cut off too low 2000 700 for RZ 1150 for PZ


% checking how many workers are currently active:
warning('off', 'MATLAB:datetime:NonstandardSystemTimeZone');
warning('off', 'MATLAB:DELETE:Permission');
p = gcp('nocreate');
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

% changing the number of active workers - if necessary:
if poolsize ~= n_parallel_workers
    if poolsize ~= 0
        delete(gcp('nocreate'));
    end
    parpool(n_parallel_workers);
end
warning('on', 'MATLAB:datetime:NonstandardSystemTimeZone');
warning('on', 'MATLAB:DELETE:Permission');
%}

% configuring opengl mode:
OS = system_dependent('getos');
if regexpi(OS, 'Windows')
    curr_opengl_mode = opengl('data');
    if curr_opengl_mode.Software && strcmp(opengl_mode, 'Hardware')
        opengl software;
    elseif ~curr_opengl_mode.Software && strcmp(opengl_mode, 'Software')
        opengl hardware;
    end
end

clc;
close all;

% asking the user to select files:
[files, directory] = get_file_name({'*.tif;*.mat'}, 'Please select the segmentation file:', 'off');
files = files{1};

if strcmp(files(end-3:end), '.mat')
    load(fullfile(directory, files));
else
    if exist(last_spacing_filename, 'file')
        load(last_spacing_filename);
    else
        last_spacing = [1 1 1];
    end
    
    % asking the user to set the spacing of the image:
    I.spacing = inputdlg( ...
        {'X spacing (microns)', 'Y spacing (microns)', 'Z spacing (microns)'}, ...
        'Please enter the spacings of the image:', ...
        1, ...
        {num2str(last_spacing(1)), num2str(last_spacing(2)), num2str(last_spacing(3))}, ...
        'on');
    
    I.spacing = cellfun(@str2num, I.spacing, 'UniformOutput', false);
    if any(cellfun(@isempty, I.spacing)) || isempty(I.spacing)
        disp('Error entering spacing values, aborting script.');
        return;
    end
    I.spacing = [I.spacing{:}];
    if any(I.spacing <= 0)
        disp('All spacing values must be positive, aborting script.');
        return;
    end
    
    last_spacing = I.spacing;
    save(last_spacing_filename, 'last_spacing');
    
    filename = fullfile(directory, files);
    I.info = imfinfo(filename);
    I.size = [I.info(1).Height, I.info(1).Width, length(I.info)];
    
    % inferring the data type:
    if I.info(1).MaxSampleValue == intmax('uint16') && I.info(1).MinSampleValue == intmin('uint16')
        I.data_type = 'uint16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('uint8') && I.info(1).MinSampleValue == intmin('uint8')
        I.data_type = 'uint8';
        I.bits = 8;
    elseif I.info(1).MaxSampleValue == intmax('int16') && I.info(1).MinSampleValue == intmin('int16')
        I.data_type = 'int16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('int8') && I.info(1).MinSampleValue == intmin('int8')
        I.data_type = 'int8';
        I.bits = 8;
    else
        disp('Data type could not be identified - contact Tomer.');
        return;
    end
    
    % allocating memory for the image field:
    I.img = zeros(I.size, I.data_type);
    
    % loading the image:
    fprintf('Loading image... ');
    for i = 1 : I.size(3)
        str = [num2str(i), '/', num2str(I.size(3))];
        fprintf(str);
        I.img(:,:,i) = imread(filename, 'Index', i);
        fprintf(repmat('\b', 1, length(str)));
    end
    disp(['Voxels = ', mat2str(I.size), ', Microns = ', mat2str(I.size .* I.spacing)]);
    
    % here we binarize the image and label each nucleus / cell, assuming 6 neighbors connectivity:
    [I.img, I.n_conncomp] = bwlabeln(I.img > 0, 6);
    disp(['Labeled binary image and identified ', num2str(I.n_conncomp), ' connected components.']);
    
    % saving the cleaned image into a mat file:
    disp('Saving image as a .mat file into the input directory..');
    save(regexprep(I.info(1).Filename, '\.[^\.]+$', '.mat'), 'I');
end

% asking the user to set the SD of the kernel and preparing the kernel:
%smooth_sd = inputdlg('Please enter the S.D. of the smoothing kernel to be used (in microns):');

%smooth_sd = str2double(smooth_sd{1});
%smooth_sd = (smooth_sd{1});
%smooth_size = round(smooth_sd ./ I.spacing * sd_to_size_coeff);
smooth_size = round(0.5 ./ I.spacing * sd_to_size_coeff);
pad_size = smooth_size + 1;
ker = fspecial3('gaussian', smooth_size);

% setting the reduction coefficient for the surfaces:
reduction_coeff = inputdlg('Please enter the coefficient of patch reduction (0, 1]:');

reduction_coeff = str2double(reduction_coeff{1});


% % asking the user whether there is an excel file with marks on the good cells:
% manual_picking = questdlg('Would you like to specify manually picked cells using an EXCEL file?', '', 'Yes', 'No', 'Yes');
% if strcmp(manual_picking, 'Yes')
%     [xls_filename, xls_dirname] = get_file_name('*.xlsx', 'Please select the excel file:', 'off');
%     I = screen_cells_using_excel_file(I, xls_dirname, xls_filename);
% end


% % here we'll identify and isolate the cells:
fprintf('Isolating cells... ');
idcs = find(I.img);
cell_idx = I.img(idcs);
[cell_idx, order] = sort(cell_idx);
idcs = idcs(order);
[x,y,z] = ind2sub(I.size, idcs);
dif = [0; find(diff(cell_idx)); length(cell_idx)];
unq_cell_idx = cell_idx(dif(1:end-1)+1);
clear idcs order;

coords = zeros(I.n_conncomp, 6, 'uint16');
bw = cell(I.n_conncomp,1);
for i = 1 : I.n_conncomp
    str = [num2str(i), '/', num2str(I.n_conncomp)];
    fprintf(str);
    rng = dif(i)+1:dif(i+1);
    coords(i,:) = [min(x(rng)), max(x(rng)), min(y(rng)), max(y(rng)), min(z(rng)), max(z(rng))];
    bw{i,1} = I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) == cell_idx(dif(i)+1);
    cc = bwconncomp(bw{i,1});
    bw{i} = imopen(bw{i}, true(3,3,3));
    fprintf(repmat('\b', 1, length(str)));
end
clear dif cell_idx x y z rng;
% I = rmfield(I, 'img');
disp(['Isolated ', num2str(length(bw)), ' cells']);

fprintf('Checking for elements with volume that is outside the specified range... ');
volume = cellfun(@nnz, bw) * prod(I.spacing);
elements_to_remove = volume < min_volume | volume > max_volume;
if nnz(elements_to_remove) > 0
    bw = bw(~elements_to_remove);
    coords(elements_to_remove,:) = [];
    unq_cell_idx = unq_cell_idx(~elements_to_remove,1);
    disp(['identified and removed ', num2str(nnz(elements_to_remove)), ' elements. Current number: ', num2str(length(bw))]);
else
    disp('none were identified.');
end
clear elements_to_remove volume;

fprintf('Removing cells that overlap with the borders of the image... ');
initial_num_cells = size(coords,1);
% non_edge_nuclei = find(~(coords(:,1) == 1 | coords(:,2) == I.size(1) | coords(:,3) == 1 | coords(:,4) == I.size(2) | coords(:,5) == 1 | coords(:,6) == I.size(3)));
% prev_unq_cell_idx = unq_cell_idx;
% unq_cell_idx = unq_cell_idx(non_edge_nuclei);
frame = false(size(I.img));
frame([1,end],:,:) = 1;
frame(:,[1,end],:) = 1;
frame(:,:,[1,end]) = 1;
is_bw_outside = false(length(bw),1);
for i = 1 : size(coords,1)
    curr_region = frame(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6));
    is_bw_outside(i) = any(curr_region(bw{i}));
    if is_bw_outside(i)
        curr_region = I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6));
        curr_region(bw{i}) = 0;
        I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) = curr_region;
    end
end
coords(is_bw_outside,:) = [];
bw(is_bw_outside) = [];
disp(['Removed ', num2str(initial_num_cells - length(bw)), ' (', num2str(round((initial_num_cells - length(bw)) / initial_num_cells * 1000) / 10), '%) cells, ' ...
    'current number of cells: ', num2str(length(bw))]);

fprintf('Extracting the surface of each cell/nucleus... ');
surfaces = repmat(struct('faces', [], 'vertices', []), length(bw), 1);
str = '';
is_problematic_nucleus = false(length(bw),1);
for i = 1 : length(bw)
    
    fprintf(repmat('\b', 1, length(str)));
    str = [num2str(i), '/', num2str(length(bw))];
    fprintf(str);
    
    padded_bw = padarray(bw{i,1}, pad_size, 0, 'both');
    smoothed = imfilter(double(padded_bw), ker);
    smoothed(smoothed < 0.001) = 0;
    
    [h, ~] = hist(nonzeros(smoothed(:)), 1000);
    cumulative_sum = cumsum(h(end:-1:1))';
    cumulative_sum = cumulative_sum(end:-1:1);
    
    [~, min_ind] = min(abs(nnz(bw{i,1}) - cumulative_sum));
    iso_value = min_ind / 1000;
    
    bw{i,1} = smoothed(pad_size(1)+1:end-pad_size(1), pad_size(2)+1:end-pad_size(2), pad_size(3)+1:end-pad_size(3)) >= iso_value;
    surface = isosurface(smoothed, iso_value);
    try
        surface.vertices = surface.vertices(:,[2,1,3]);
    catch
        is_problematic_nucleus(i) = true;
        continue;
    end
    
    surface.vertices = (surface.vertices - repmat(pad_size, size(surface.vertices,1), 1) + repmat(double(coords(i,[1,3,5]))-1, size(surface.vertices,1), 1)) .* repmat(I.spacing, size(surface.vertices,1), 1); 
    
    surface = reducepatch(surface, reduction_coeff);
    
    surfaces(i,1) = surface;
    
    if 0
        figure('Color', [0 0 0]); %#ok<*UNRCH>
        patch(surface, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        rotate3d on;
        daspect([1,1,1])
        axis vis3d;
        camlight;
        axis off;
    end
end

n_problematic_nuclei = nnz(is_problematic_nucleus);
if n_problematic_nuclei > 0
    coords(is_problematic_nucleus,:) = [];
    bw(is_problematic_nucleus) = [];
    surfaces(is_problematic_nucleus) = [];
    disp(['Removed ', num2str(n_problematic_nuclei), ' nuclei, total remaining nuclei: ', num2str(length(bw))]);
end

% saving the surfaces into a .mat file as a structure array:
spacing = I.spacing; %#ok<*NASGU>
save(fullfile(directory, regexprep(files, '\..*$', ' (Surfaces).mat')), 'surfaces', 'bw', 'coords', 'spacing');

% joining all cells/nuclei into the same image:
vertices_len = arrayfun(@(x) length(x.vertices), surfaces);
faces_len = arrayfun(@(x) length(x.faces), surfaces);
vertices_cumlen = [0; cumsum(vertices_len(1:end))];
faces_cumlen = [0; cumsum(faces_len(1:end))];

% updating the indices of the faces in each induvidual surface:
for i = 2 : length(surfaces)
    surfaces(i).faces = surfaces(i).faces + vertices_cumlen(i);
end

% allocating memory for the joint surface:
S.vertices = nan(sum(vertices_len),3);
S.faces = nan(sum(faces_len),3);

% joining the surfaces:
for i = 1 : length(vertices_cumlen)-1
    S.vertices(vertices_cumlen(i)+1:vertices_cumlen(i+1),:) = surfaces(i).vertices;
    S.faces(faces_cumlen(i)+1:faces_cumlen(i+1),:) = surfaces(i).faces;
end

% saving the joined surface:
if 0
    fprintf('Saving the output file of the surfaces... ');
    S.vertices = S.vertices / 1000;
    % S.vertices = S.vertices(:,[2,1,3]);
    write_vtk(fullfile(directory, regexprep(files, '\..*$', ' (Surfaces).vtk')), S.vertices, S.faces);
    disp('done.');
end

% creating a .vff image for the screened segmented cells and saving it:
if 0
    fprintf('Saving the output file of the volume... ');
    J.header = vffread(demo_vff_file, '', 1, 'header');
    required_fields = {'format', 'miscellaneous', 'rank', 'scanner', 'Z1', 'type'};
    fn = fieldnames(J.header);
    J.header = rmfield(J.header, setdiff(fn, required_fields));
    J.img = I.img;
    J.size = size(J.img);
    J.header.size = size(J.img);
    J.header.rawsize = numel(J.img) * I.bits / 8;
    J.data_type = I.data_type;
    J.header.data_type = I.data_type;
    J.header.bits = I.bits;
    J.header.origin_in_mm = [0, 0, 0];
    J.header.origin_in_pix = [0, 0, 0];
    J.spacing = I.spacing / 1000;
    J.header.spacing = I.spacing / 1000;
    output_filename = fullfile(directory, regexprep(files, '\..*$', ' (Screened).vff'));
    vffwrite(J, output_filename);
    disp('done.');
end

show_surface(S);
cameramenu;
light;

if 0
    figure('Color', [0 0 0]);
    patch(S, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    rotate3d on;
    daspect([1,1,1])
    axis vis3d;
    camlight;
    axis off;
end

fprintf('\nDone.\n');

return;

[vol, bw, unq_cell_idx, coords] = calc_nuclei_volumes(bw, I.spacing, unq_cell_idx, coords);
[vol, coords, bw, unq_cell_idx] = cutoff_volume(vol, coords, bw, unq_cell_idx);

%%%%

% fprintf('Calculating the surface of each cell... ');
% disp('Done!');

fprintf('Calculating the volume and surface area of each cell... ');
surface_volume = zeros(length(surfaces), 1, 'single');
surface_area = zeros(length(surfaces), 1, 'single');
for i = 1 : length(surfaces)
    str = [num2str(i), '/', num2str(length(surfaces)-1)];
    fprintf(str);
    e1 = surfaces{i}.vertices(surfaces{i}.faces(:,1),:) - surfaces{i}.vertices(surfaces{i}.faces(:,2),:);
    e3 = surfaces{i}.vertices(surfaces{i}.faces(:,3),:) - surfaces{i}.vertices(surfaces{i}.faces(:,1),:);
    [surface_area(i,1),surface_volume(i,1)] = stlVolumeNormals(surfaces{i}.vertices', surfaces{i}.faces', cross(e1,e3));
    fprintf(repmat('\b', 1, length(str)));
end
disp('Done!');

keyboard;

% figure
% subplot(1,2,1)
% fv = isosurface(data,.5);
% p1 = patch(fv,'FaceColor','red','EdgeColor','none');
% view(3)
% daspect([1,1,1])
% axis tight
% camlight
% camlight(-80,-10)
% lighting gouraud
% title('Triangle Normals')
%
% subplot(1,2,2)
% fv = isosurface(data,.5);
% p2 = patch(fv,'FaceColor','red','EdgeColor','none');
% n = isonormals(data,p2);
% view(3)
% daspect([1 1 1])
% axis tight
% camlight
% camlight(-80,-10)
% lighting gouraud
% title('Data Normals')

%%%%

% saving volumes:
for i = 1 : length(unq_cell_idx)
    J.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) = uint16(bw{i,1}) * uint16(round(vol(i)));
end
vffwrite(J, output_filename);


function [vol, bw, unq_cell_idx, coords] = calc_nuclei_volumes(bw, spacing, unq_cell_idx, coords)

% calculating the volume of each cell and sorting the cells based on their volumes:
vol = cellfun(@nnz, bw) * prod(spacing);
[vol, order] = sort(vol);
unq_cell_idx = unq_cell_idx(order);
coords = coords(order,:);
bw = bw(order);


function [vol, coords, bw, unq_cell_idx] = cutoff_volume(vol, coords, bw, unq_cell_idx) %#ok<*DEFNU>

% plotting the volumes and the distribution of the volumes:
figure('Position', [279 277 1231 598]);
subplot(1,2,1);
plot(vol, '.');
xlabel('Nucleus serial number');
ylabel('Nucleus volume (micron^3)');
title('Nuclei volume');
axis square;
subplot(1,2,2);
hist(vol, sqrt(length(vol)));
xlabel('Nucleus volume (micron^3)');
ylabel('# Appearances');
title('Nuclei volume distribution');
axis square;

set_volume_cutoff = questdlg('Would you like to set a cutoff for the volume of the nuclei?', '', 'Yes - based on volumes', 'Yes - based on volume distribution', 'No', 'Yes - based on volumes');

% processing user response:
if isempty(set_volume_cutoff)
    error('No value was selected, aborting script.');
end

ax = get(gcf, 'Children');
if strcmp(set_volume_cutoff, 'Yes - based on volumes')
    axes(ax(2));
    vol_range = getrect;
    vol_range = vol_range([1,3]);
    idcs_to_take = ceil(vol_range(1)) : floor(vol_range(2));
elseif strcmp(set_volume_cutoff, 'Yes - based on volume distribution')
    axes(ax(1));
    vol_range = getrect;
    vol_range = vol_range([1,3]);
    idcs_to_take = find(vol >= vol_range(1) & vol <= vol_range(2));
else
    idcs_to_take = 1:length(vol);
end

initial_num_cells = length(vol);
vol = vol(idcs_to_take);
coords = coords(idcs_to_take,:);
bw = bw(idcs_to_take);
unq_cell_idx = unq_cell_idx(idcs_to_take);
disp(['Removed ', num2str(initial_num_cells - length(bw)), ' (', num2str(round((initial_num_cells - length(bw)) / initial_num_cells * 1000) / 10), '%) cells, ' ...
    'current number of cells: ', num2str(length(bw))]);


function I = screen_cells_using_excel_file(I, xls_dirname, xls_filename)

[P, ~, C] = xlsread(fullfile(xls_dirname, xls_filename{1}));
x_col = find(strcmp(C(1,:), 'X'));
y_col = find(strcmp(C(1,:), 'Y'));
z_col = find(strcmp(C(1,:), 'Slice'));

p = round(P(:,[x_col, y_col, z_col]));
to_remove = p(:,1) > I.size(2) | p(:,2) > I.size(1) | p(:,3) > I.size(3) | any(p < 1, 2);
p(to_remove,:) = [];

I.cell_idx = unique(nonzeros(I.img(sub2ind(I.size, p(:,2), p(:,1), p(:,3)))));
I.img(~ismember(I.img, I.cell_idx)) = 0;


function I = find_repeating_indices(I)

% here we'll verify that the indexes of the cells are unique:
fprintf('Identifying cells with unique indices... ');
idcs = find(I.img);
cell_idx = I.img(idcs);
[cell_idx, order] = sort(cell_idx);
idcs = idcs(order);
[x,y,z] = ind2sub(I.size, idcs);
dif = [0; find(diff(cell_idx)); length(cell_idx)];
unq_cell_idx = cell_idx(dif(1:end-1)+1);
disp(['Identified ', num2str(length(unq_cell_idx)), ' unique cells']);
clear idcs order;

coords = zeros(length(dif)-1, 6, 'uint16');
non_unique_idcs = nan(2,0);
last_idx = unq_cell_idx(end);
bw = cell(length(dif)-1,1);
for i = 1 : length(dif)-1
    %     str = [num2str(i), '/', num2str(length(dif)-1)];
    %     fprintf(str);
    rng = dif(i)+1:dif(i+1);
    coords(i,:) = [min(x(rng)), max(x(rng)), min(y(rng)), max(y(rng)), min(z(rng)), max(z(rng))];
    bw{i,1} = I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) == cell_idx(dif(i)+1);
    cc = bwconncomp(bw{i,1});
    if cc.NumObjects > 1
        cc = bwconncomp(I.img == cell_idx(dif(i)+1));
        non_unique_idcs(:,end+1:end+cc.NumObjects) = [repmat(cell_idx(dif(i)+1), 1, cc.NumObjects); cell_idx(dif(i)+1), last_idx+1:last_idx+cc.NumObjects-1];
        for j = 2 : cc.NumObjects
            last_idx = last_idx + 1;
            I.img(cc.PixelIdxList{j}) = last_idx;
        end
    end
    %     fprintf(repmat('\b', 1, length(str)));
end

I.non_unique_idcs = non_unique_idcs;

if ~isempty(non_unique_idcs)
    disp(['NOTE: Identified ', num2str(size(non_unique_idcs,2)), ' cells/nuclei with non-unique indices, ', ...
        'saving corrected image as a .mat file in the same directory:']);
    disp(non_unique_idcs);
end

% if ~isempty(non_unique_idcs)
%     questdlg(['NOTE: Identified ', num2str(size(non_unique_idcs,2)), ' cells/nuclei with non-unique indices (see Command Window for details), ', ...
%         'would you like to save the corrected image?']);
%     [FileName, PathName, FilterIndex] = uiputfile(',DialogTitle);
%
% end

