function plot_cell_morphology_analysis

%{
INSTRUCTIONS:
This script receives from the user the four .mat files that contain the characteristics of the four growth plate zones, calculates some additional
characteristics, and finally plots all the graphs for the user.
%}

% closing all open figures:
close all;
ignore_cells = true;
wait_after_view_link = false;
opengl_mode = 'Hardware'; % opengl_mode can be either 'Hardware' (somewhat slower but complex figures might be in higher quality), or 'Software' (faster, but complex
% figures might be in lower quality). The selection of the mode of opengl can be done this way only in computers that run Windows as the OS, for other OS please
% see documentation for opengl (write 'doc opengl' in the command window).

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

% ask user to select the characteristics files:
[filename, dirname] = get_file_name('*.mat', 'Please select all the characteristics files:', 'on');

% loading all the files:
for i = 1 : length(filename)
    zone_name = regexprep(filename{i}, ' \(Characteristics\).mat', '');
    data = load(fullfile(dirname, filename{i}));
    G.(zone_name) = data.G; % general
    if ~ignore_cells
        C.(zone_name) = data.C; % cells
    end
    N.(zone_name) = data.N; % nuclei
end

% adding surface area to volume ratio to the calculated characteristics:
G = add_characteristics(G, ignore_cells);

% here we start plotting all the resutls:
% if ~ignore_cells
%     plot_pca(G, 'cel');
% end
% plot_pca(G, 'nuc');

if ~ignore_cells
%     plot_ellipsoid_approximations(G, 'cel', 'ellipsoid_radii');
    plot_ellipsoid_approximations(G, 'cel', 'PC_range');
end
% plot_ellipsoid_approximations(G, 'nuc', 'ellipsoid_radii');
plot_ellipsoid_approximations(G, 'nuc', 'PC_range');

plot_individual_distributions(G, ignore_cells);
plot_triplet_characteristics(G, {'volume', 'surface_area', 'sphericity'}, ignore_cells);

if ~ignore_cells
    % for pc = 1 : 3
    %     plot_shape_orientations(G, 'cel', pc);
    % end
    plot_shape_orientations_arrows(G, 'cel', wait_after_view_link);
end

% for pc = 1 : 3
%     plot_shape_orientations(G, 'nuc', pc);
% end
plot_shape_orientations_arrows(G, 'nuc', wait_after_view_link);

% plot_nucleus_shift_orientations(G); % this function has not been quality assured yet...


function G = add_characteristics(G, ignore_cells)

zones = fieldnames(G);

% surface area to volume ratio:
for i = 1 : length(zones)
    if ~ignore_cells
        G.(zones{i}).nuc.sa2vol_ratio = G.(zones{i}).nuc.surface_area ./ G.(zones{i}).nuc.volume;
    end
    G.(zones{i}).nuc.sa2vol_ratio = G.(zones{i}).nuc.surface_area ./ G.(zones{i}).nuc.volume;
end


function plot_ellipsoid_approximations(G, element, radii)

zones = fieldnames(G);
colors = jet(length(zones));
colors(3,:) = [0 1 0];

vals = [];
for z = 1 : length(zones)
    vals = [vals; G.(zones{z}).(element).(radii), z*ones(size(G.(zones{z}).(element).(radii),1),1)];
end
% xl = [min(vals(:,1)), max(vals(:,1))];
% yl = [min(vals(:,2)), max(vals(:,2))];
% zl = [min(vals(:,3)), max(vals(:,3))];
% xl = [xl(1)-0.1*range(xl), xl(2)+0.1*range(xl)];
% yl = [yl(1)-0.1*range(yl), yl(2)+0.1*range(yl)];
% zl = [zl(1)-0.1*range(zl), zl(2)+0.1*range(zl)];

% first, plotting pca with original scale:
views = [0, 90; 0, 0; 90, 0; 19, 25];

% -----------------------------------------------------------------------------------------

fid = figure('Name', [element, ' - Radii of fitted ellipsoid']);
for i = 1 : 4
    subplot(2,2,i)
    hold on;
    for z = 1 : length(zones)
        idcs = vals(:,4) == z;
        scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
    end
    hold off;
    axis image;
    grid on;
    xlabel('Axis 1');
    ylabel('Axis 2');
    zlabel('Axis 3');
    view(views(i,1), views(i,2));
    
    legend(zones, 'Location', 'BestOutside');
end
xl = xlim; xlim([0, xl(2) + range(xl)*0.1]);
yl = ylim; ylim([0, yl(2) + range(yl)*0.1]);
zl = zlim; zlim([0, zl(2) + range(zl)*0.1]);
rotate3d on;
axis vis3d;

print('-RGBImage');

% -----------------------------------------------------------------------------------------

original_vals = vals;
vals(:,1:3) = vals(:,1:3) ./ repmat(sum(vals(:,1:3),2),1,3);

fid = figure('Name', [element, ' - Normalized radii of fitted ellipsoid']);
hold on;
for z = 1 : length(zones)
    idcs = vals(:,4) == z;
    scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
end
hold off;
grid on;
xlabel('Axis 1');
ylabel('Axis 2');
zlabel('Axis 3');
axis image;
view(125,35);

xl = xlim; xlim([xl(1) - range(xl)*0.1, xl(2) + range(xl)*0.1]);
yl = ylim; ylim([yl(1) - range(yl)*0.1, yl(2) + range(yl)*0.1]);
zl = zlim; zlim([zl(1) - range(zl)*0.1, zl(2) + range(zl)*0.1]);
rotate3d on;
set(gca, 'XDir', 'reverse');

print('-RGBImage');

% -----------------------------------------------------------------------------------------

movingPoints = [1 0; 0.5 0.5; 1/3 1/3];
fixedPoints = [0 0; 0 1; 1 1];
tform = fitgeotrans(movingPoints, fixedPoints, 'Affine');
vals = [vals(:,1:2), ones(size(vals,1),1)] * tform.T;
vals(:,3) = [];

% vals(:,2) = vals(:,2) / 3;

fid = figure('Name', [element, ' - Aspect ratio profile']);
subplot(2,4,5);
hold on;
for z = 1 : length(zones)
    idcs = original_vals(:,4) == z;
    scatter(vals(idcs,2), vals(idcs,1), 5, 'filled', 'MarkerFaceColor', colors(z,:));
end
plot([0,1], [0,1], 'Color', [1,1,1]*0.85);
hold off;
axis image;
xlim([0,1]);
ylim([0,1]);
% ylabel('PC1 % Variability');
% xlabel('PC2 % Variability');
set(gca, 'YAxisLocation', 'right');
legend(zones, 'Location', 'NorthWest');

for z = 1 : length(zones)
    subplot(2,4,z);
    idcs = original_vals(:,4) == z;
    scatter(vals(idcs,2), vals(idcs,1), 5, 'filled', 'MarkerFaceColor', colors(z,:));
    hold on;
    plot([0,1], [0,1], 'Color', [1,1,1]*0.85);
    hold off;
    axis image;
    xlim([0,1]);
    ylim([0,1]);
    legend(zones{z}, 'Location', 'NorthWest');
end
% ylabel('PC1 % Variability');
% xlabel('PC2 % Variability');
set(gca, 'YAxisLocation', 'right');

% here we'll calculate the convex hull of the points of each zone and plot them together:
subplot(2,4,6);
hold on;
for z = length(zones):-1:1
    idcs = find(original_vals(:,4) == z);
    K = convhull(vals(idcs,2), vals(idcs,1));
    plot(vals(idcs(K),2), vals(idcs(K),1), '--', 'Color', colors(z,:), 'LineWidth', 1);
end
plot([0,1], [0,1], 'Color', [1,1,1]*0.85);
hold off;
axis image;
xlim([0,1]);
ylim([0,1]);
set(gca, 'YAxisLocation', 'right');

print('-RGBImage');

% archetype points:
% cells top right: #391, 176 in rz
% show_surface(data.rz.surfaces(176));
% cells bottom left: #190, 77 in pz
% cells bottom right: #37, 20 in phz
% cells center: #168, 55 in pz


function plot_pca(G, entity)

zones = fieldnames(G);
colors = jet(length(zones));
colors(3,:) = [0 1 0];

vals = [];
for z = 1 : length(zones)
    vals = [vals; G.(zones{z}).(entity).ellipsoid_radii, z*ones(size(G.(zones{z}).(entity).ellipsoid_radii,1),1)];
end
% xl = [min(vals(:,1)), max(vals(:,1))];
% yl = [min(vals(:,2)), max(vals(:,2))];
% zl = [min(vals(:,3)), max(vals(:,3))];
% xl = [xl(1)-0.1*range(xl), xl(2)+0.1*range(xl)];
% yl = [yl(1)-0.1*range(yl), yl(2)+0.1*range(yl)];
% zl = [zl(1)-0.1*range(zl), zl(2)+0.1*range(zl)];

% first, plotting pca with original scale:
views = [0, 90; 0, 0; 90, 0; 19, 25];

% ---------------------------------------------------------------------------------------------------

if strcmp(entity, 'cel')
    figure('Name', 'Cells - absolute size of axes');
else
    figure('Name', 'Nuclei - absolute size of axes');
end

for i = 1 : 4
    subplot(2,2,i)
    hold on;
    for z = 1 : length(zones)
        idcs = vals(:,4) == z;
        scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
    end
    hold off;
    axis image;
    axis vis3d;
    grid on;
    xlabel('Axis 1');
    ylabel('Axis 2');
    zlabel('Axis 3');
    view(views(i,1), views(i,2));
    
    xl = xlim; xlim([0, xl(2)*1.15]);
    yl = ylim; ylim([0, yl(2)*1.15]);
    zl = zlim; zlim([0, zl(2)*1.15]);
end
rotate3d on;
axis vis3d;

% ---------------------------------------------------------------------------------------------------

original_vals = vals;
vals(:,1:3) = vals(:,1:3) ./ repmat(sum(vals(:,1:3),2),1,3);

if strcmp(entity, 'cel')
    figure('Name', 'Cells - relative size of axes');
else
    figure('Name', 'Nuclei - relative size of axes');
end

hold on;
for z = 1 : length(zones)
    idcs = vals(:,4) == z;
    scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
end
hold off;
grid on;
xlabel('Axis 1');
ylabel('Axis 2');
zlabel('Axis 3');
axis image;
view(125,35);
xl = xlim; xlim([xl(1) - range(xl)*0.1, xl(2) + range(xl)*0.1]);
yl = ylim; ylim([yl(1) - range(yl)*0.1, yl(2) + range(yl)*0.1]);
zl = zlim; zlim([zl(1) - range(zl)*0.1, zl(2) + range(zl)*0.1]);
rotate3d on;
axis vis3d;
set(gca, 'XDir', 'reverse');

% ---------------------------------------------------------------------------------------------------

movingPoints = [1 0; 0.5 0.5; 1/3 1/3];
fixedPoints = [0 0; 0 1; 1 1];
tform = fitgeotrans(movingPoints, fixedPoints, 'Affine');
vals = [vals(:,1:2), ones(size(vals,1),1)] * tform.T;
vals(:,3) = [];

% vals(:,2) = vals(:,2) / 3;

if strcmp(entity, 'cel')
    figure('Name', 'Cells dimensions');
else
    figure('Name', 'Nuclei dimensions');
end

hold on;
for z = 1 : length(zones)
    idcs = original_vals(:,4) == z;
    scatter(vals(idcs,2), vals(idcs,1), 5, 'filled', 'MarkerFaceColor', colors(z,:));
end
plot([0,1], [0,1], 'Color', [1,1,1]*0.85);
hold off;
axis image;
xlim([0,1]);
ylim([0,1]);
% ylabel('PC1 % Variability');
% xlabel('PC2 % Variability');
set(gca, 'YAxisLocation', 'right');
legend(regexprep(zones, '_', ' '), 'Location', 'NorthWest');


function plot_nucleus_shift_orientations(G)

zones = fieldnames(G);

vals = [];
for z = 1 : length(zones)
    vals = [vals; G.(zones{z}).inter.centroid_shift_physical_coords];
end
xl = [min(vals(:,1)), max(vals(:,1))];
yl = [min(vals(:,2)), max(vals(:,2))];
zl = [min(vals(:,3)), max(vals(:,3))];
xl = [xl(1)-0.1*range(xl), xl(2)+0.1*range(xl)];
yl = [yl(1)-0.1*range(yl), yl(2)+0.1*range(yl)];
zl = [zl(1)-0.1*range(zl), zl(2)+0.1*range(zl)];

figure;
colors = jet(length(zones));
colors(3,:) = [0 0 0];

for z = 1 : length(zones)
    hax(z) = subplot(2,3,z);
    vals = [G.(zones{z}).inter.centroid_shift_physical_coords, z * ones(G.(zones{z}).inter.n_overlapping,1)];
    %     scatter3(vals(:,1), vals(:,2), vals(:,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
    quiver3(zeros(size(vals,1),1), zeros(size(vals,1),1), zeros(size(vals,1),1), vals(:,1), vals(:,2), vals(:,3), 1, 'Color', colors(z,:));
    axis image;
    axis vis3d;
    rotate3d on;
    xlim(xl);
    ylim(yl);
    zlim(zl);
    grid off;
end

hax(5) = subplot(2,3,5);
hold on;
for z = 1 : length(zones)
    vals = [G.(zones{z}).inter.centroid_shift_physical_coords, z * ones(G.(zones{z}).inter.n_overlapping,1)];
    %     scatter3(vals(:,1), vals(:,2), vals(:,3), 5, 'filled', 'MarkerFaceColor', colors(z,:));
    quiver3(zeros(size(vals,1),1), zeros(size(vals,1),1), zeros(size(vals,1),1), vals(:,1), vals(:,2), vals(:,3), 1, 'Color', colors(z,:));
end
hold off;
xlim(xl);
ylim(yl);
zlim(zl);
axis image;
axis vis3d;
rotate3d on;
grid off;
legend(zones, 'Location', 'NorthEastOutside');

linkprop(hax, 'view');


function plot_shape_orientations_arrows(G, entity, wait_after_view_link)

% preparing parameters:
zones = fieldnames(G);
colors = [1 0 0; 0 1 0; 0 0 1];

% collecting the orientation vectors of all cells / nuclei in each zone:
vals = cell(1,3);
for pc = 1 : 3
    for z = 1 : length(zones)
        tmp = cell2mat(cellfun(@(x) x(:,pc)', G.(zones{z}).(entity).PCA_coeff, 'UniformOutput', false));
        tmp(:,4) = z;
%         tmp(:,5) = sum(G.(zones{z}).(entity).ellipsoid_radii(:,1:2),2) ./ sum(G.(zones{z}).(entity).ellipsoid_radii,2);
        vals{pc} = [vals{pc}; tmp];
    end
end

% making sure that the selected directions of all orientation vectors (out of the two) are at the same hemisphere:
for pc = [1, 3]
    invert_orientation = vals{pc}(:,4-pc) < 0;
    vals{pc}(invert_orientation, 1:3) = -vals{pc}(invert_orientation, 1:3);
end
vals{2}(:,1:3) = cross(vals{1}(:,1:3), vals{3}(:,1:3));

% Plotting one sphere per zone with the orientation of its cells / nuclei:
arrow_start_coeff = 0.8;
axis_half_length = 1.2;
figure('Color', 'w', 'Name', 'Orientations of PC1 (red) and PC3 (blue)');
for z = 1 : length(zones)
    hax(z) = subplot(1,length(zones),z); %#ok<AGROW>
    hold on;
    for pc = [1,3]
        idcs = vals{pc}(:,4) == z;
        
        curr_vals = nan(nnz(idcs)*3-1,3);
        curr_vals(1:3:end,:) = vals{pc}(idcs,1:3) * arrow_start_coeff;
        curr_vals(2:3:end,:) = vals{pc}(idcs,1:3);
        
        plot3(curr_vals(:,1), curr_vals(:,2), curr_vals(:,3), 'Color', colors(pc,:));
        
        if pc == 1
            quiver3([-axis_half_length;0;0], [0;-axis_half_length;0], [0;0;-axis_half_length], ...
                [axis_half_length;0;0], [0;axis_half_length;0], [0;0;axis_half_length], 2, 'Color', [1 1 1]*0.75);
            
            text(axis_half_length * 1.15, 0, 0, '+X', 'HorizontalAlignment', 'center');
            text(-axis_half_length * 1.15, 0, 0, '-X', 'HorizontalAlignment', 'center');
            text(0, axis_half_length * 1.15, 0, '+Y', 'HorizontalAlignment', 'center');
            text(0, -axis_half_length * 1.15, 0, '-Y', 'HorizontalAlignment', 'center');
            text(0, 0, axis_half_length * 1.15, '+Z', 'HorizontalAlignment', 'center');
            text(0, 0, -axis_half_length * 1.15, '-Z', 'HorizontalAlignment', 'center');
            axis image;
            rotate3d on;
            grid off;
            axis off;
        end
    end
    hold off;
    title(regexprep(zones{z}, '_', ' '));
end

for i = 1 : length(zones)
    subplot(1,length(zones),i);
    axis vis3d;
end

cameratoolbar;

print('-RGBImage');

if wait_after_view_link
    linkprop(hax, 'view');
    disp('To continue press the F5 key.');
    keyboard;
end

return;

% here we calculate the spheric density of each principal component:
R = 1; % the radius of the sphere
ST_area = cell(3,1);
point_area = cell(3,1);
point_density = cell(3,1);
m = cell(3,1);
for pc = 1 : 3
    
    % first, since there may be more than one cell or nucleus to the same direction, we assign a weight to each cell / nucleus:
    unq_vals = unique(vals{pc}(:,1:3), 'rows');
    m{pc} = zeros(size(unq_vals,1),1);
    for i = 1 : size(unq_vals,1)
        m{pc}(i,1) = nnz(vals{pc}(:,1) == unq_vals(i,1) & vals{pc}(:,2) == unq_vals(i,2) & vals{pc}(:,3) == unq_vals(i,3));
    end
    
    curr_vals = [unq_vals; -unq_vals];
    m{pc} = [m{pc}; m{pc}];
    
    [~, face] = sphere_delaunay(size(curr_vals,1), curr_vals');
    face = face';
    
    % checking if there are triangles that were not assigned a face:
    % first, calculating how many neighbors each vertex has:
    neigh = cell(size(curr_vals,1),1);
    n_triangles = zeros(size(curr_vals,1),1);
    for i = 1 : size(curr_vals,1)
        neigh{i} = setdiff(unique(face(any(face == i,2),:)), i);
        n_triangles(i,1) = nnz(any(face == i,2));
    end
    if ~all(cellfun(@length, neigh) == n_triangles)
        disp('Unreferenced triangles found, contact Tomer.');
        keyboard;
    end
    
    ST_area{pc} = zeros(size(face,1),1);
    
    for f = 1 : size(face,1)
        
        % the three points of the triangle:
        pt1 = curr_vals(face(f,1),:);
        pt2 = curr_vals(face(f,2),:);
        pt3 = curr_vals(face(f,3),:);
        
        % the spherical lengths of the three edges of the triangle:
        a = acos(dot(pt1, pt2));
        b = acos(dot(pt1, pt3));
        c = acos(dot(pt2, pt3));
        
        % the semiperimeter:
        s = (a + b + c) / 2;
        
        % the spherical excess:
        E = 4 * atan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)));
        
        % the area of the triangle:
        ST_area{pc}(f,1) = R^2 * E;
        
    end
    
    for i = 1 : size(curr_vals,1)
        point_area{pc}(i,1) = sum(ST_area{pc}(any(face == i,2)));
    end
    
    point_density{pc} = 3 * m{pc} ./ point_area{pc};
    
    % now we add more interpolated points:
    n_rand = round(ST_area{pc} ./ prctile(ST_area{pc}, 5));
    for f = 1 : size(face,1)
        pt1 = curr_vals(face(f,1),:);
        pt2 = curr_vals(face(f,2),:);
        pt3 = curr_vals(face(f,3),:);
        
        r = rand(n_rand(f), 2);
        switch_idcs = r(:,1) < r(:,2);
        r(switch_idcs,:) = 1-r(switch_idcs,:);
        r_points = repmat(pt1, size(r,1), 1) + ...
            repmat(r(:,1), 1, 3) .* repmat(pt2-pt1, size(r,1), 1) + ...
            repmat(r(:,2), 1, 3) .* repmat(pt3-pt2, size(r,1), 1);
        
        % interpolating the colors to each random point:
        for i = 1 : size(r_points,1)
            T1_area = sqrt(det([r_points(i,:); pt2; ones(1,3)])^2 + det([pt2; pt3; ones(1,3)])^2 + det([pt3; r_points(i,:); ones(1,3)])^2) / 2;
            T2_area = sqrt(det([r_points(i,:); pt1; ones(1,3)])^2 + det([pt1; pt3; ones(1,3)])^2 + det([pt3; r_points(i,:); ones(1,3)])^2) / 2;
            T3_area = sqrt(det([r_points(i,:); pt1; ones(1,3)])^2 + det([pt1; pt2; ones(1,3)])^2 + det([pt2; r_points(i,:); ones(1,3)])^2) / 2;
        end
        
        keyboard;
        Vq = interp3(X, Y, Z, point_density{1}(face(f,:)), r_points(:,1), r_points(:,2), r_points(:,3));
    end
    
    [x,y,z] = sphere(1000);
    [s_azimuth, s_elevation, s_r] = cart2sph(x,y,z);
    clear x y z;
    [v_azimuth, v_elevation, v_r] = cart2sph(curr_vals(:,1), curr_vals(:,2), curr_vals(:,3));
    
    all_v_azimuth = v_azimuth;
    all_v_elevation = v_elevation;
    
    all_v_azimuth = [all_v_azimuth; -v_azimuth];
    all_v_elevation = [all_v_elevation; pi-v_elevation];
    
    F = scatteredInterpolant(v_azimuth, v_elevation, log(point_density{pc}));
    interp_density{pc} = F(s_azimuth, s_elevation);
    
    % transferring back to the original coordinate system:
    [x,y,z] = sph2cart(s_azimuth, s_elevation, s_r);
    clear s_azimuth s_elevation s_r;
    
    close all;
    figure;
    surf(x, y, z, interp_density{pc}, 'EdgeColor', 'none', 'FaceColor', 'interp');
    colormap jet;
    colorbar;
    axis image;
    axis vis3d;
    rotate3d on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end


function plot_shape_orientations(G, entity, pc)

% preparing parameters:
zones = fieldnames(G);
colors = jet(length(zones));
colors(3,:) = [0 0 0];
[X,Y,Z] = sphere(100);
markers = {'s', '+', '.'};
marker_size = [20, 20, 30];

figure('Name', [entity, ' PC', num2str(pc)]);
for z = 1 : length(zones)+1
    subplot(1, length(zones)+1, z);
    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', [1,1,1]*0.5, 'FaceAlpha', 0.1);
end

% collecting the orientation vectors of all cells / nuclei in each zone:
vals = [];
for z = 1 : length(zones)
    tmp = cell2mat(cellfun(@(x) x(:,pc)', G.(zones{z}).(entity).PCA_coeff, 'UniformOutput', false));
    tmp(:,4) = z;
    tmp(:,5) = sum(G.(zones{z}).(entity).ellipsoid_radii(:,1:2),2) ./ sum(G.(zones{z}).(entity).ellipsoid_radii,2);
    vals = [vals; tmp]; %#ok<AGROW>
end

% making sure that the selected directions of all orientation vectors (out of the two) are at the same hemisphere:
invert_orientation = vals(:,4-pc) < 0;
vals(invert_orientation, 1:3) = -vals(invert_orientation, 1:3);

% adding to the same figure one sphere per zone with the orientation of its cells / nuclei:
for z = 1 : length(zones)
    hax(z) = subplot(1,length(zones)+1,z); %#ok<AGROW>
    idcs = vals(:,4) == z;
    hold on;
    scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), ...
        markers{pc}, 'SizeData', marker_size(pc), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors(z,:));
    hold off;
    axis image;
    axis vis3d;
    xlim([-1,1]);
    ylim([-1,1]);
    zlim([-1,1]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    rotate3d on;
    grid off;
end

% plotting a sphere with the orientations of all zones together:
hax(5) = subplot(1,length(zones)+1,5);
hold on;
for z = 1 : length(zones)
    idcs = vals(:,4) == z;
    scatter3(vals(idcs,1), vals(idcs,2), vals(idcs,3), ...
        markers{pc}, 'SizeData', marker_size(pc), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors(z,:));
end
hold off;
axis image;
axis vis3d;
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
xlabel('X');
ylabel('Y');
zlabel('Z');
rotate3d on;
grid off;
% legend(zones, 'Location', 'NorthEastOutside');

% linking the orientations of all spheres together and allowing the user to keyboard and select the best view:
linkprop(hax, 'view');
% keyboard;

% for a single selected zone, we plot the semi-transparent surfaces of its entities and the corresponding orientation vectors:
if 0
    load(['E:\Unsynced data\Data for Hagai\', zones{1}, '\mem_zcrop_segmented (Surfaces).mat']);
    show_surface(surfaces);
    directions = cell2mat(cellfun(@(x) x(:,pc)', G.(zones{1}).(entity).PCA_coeff, 'UniformOutput', false));
    centroids = G.(zones{1}).(entity).centroids;
    hold on;
    quiver3(centroids(:,1), centroids(:,2), centroids(:,3), directions(:,1), directions(:,2), directions(:,3), 0.5);
    hold off;
    cameramenu;
end


function plot_triplet_characteristics(G, characteristics, ignore_cells)

zones = fieldnames(G);
figure;

if ~ignore_cells
    vals = [];
    for z = 1 : length(zones)
        tmp = [G.(zones{z}).cel.(characteristics{1}), G.(zones{z}).cel.(characteristics{2}), G.(zones{z}).cel.(characteristics{3})];
        tmp(:,4) = z;
        vals = [vals; tmp];
    end
    
    hax(1) = subplot(1,2,1);
    scatter3(vals(:,1), vals(:,2), vals(:,3), 5, vals(:,4), 'filled');
    axis vis3d;
    cameratoolbar;
    xlabel(regexprep(characteristics{1}, '_', ' '));
    ylabel(regexprep(characteristics{2}, '_', ' '));
    zlabel(regexprep(characteristics{3}, '_', ' '));
    title('Cells');
end

vals = [];
for z = 1 : length(zones)
    tmp = [G.(zones{z}).nuc.(characteristics{1}), G.(zones{z}).nuc.(characteristics{2}), G.(zones{z}).nuc.(characteristics{3})];
    tmp(:,4) = z;
    vals = [vals; tmp];
end

hax(2) = subplot(1,2,2);
scatter3(vals(:,1), vals(:,2), vals(:,3), 5, vals(:,4), 'filled');
axis vis3d;
cameratoolbar;
xlabel(regexprep(characteristics{1}, '_', ' '));
ylabel(regexprep(characteristics{2}, '_', ' '));
zlabel(regexprep(characteristics{3}, '_', ' '));
title('Nuclei');

cameratoolbar;

print('-RGBImage');

if ~ignore_cells && wait_after_view_link
    linkprop(hax, 'view');
    disp('To continue press the F5 key.');
    keyboard;
end


function plot_individual_distributions(G, ignore_cells)

% plotting the distributions of the different characteristics:
fid = figure;
idx = 1;

zones = fieldnames(G);
characteristics = setdiff(fieldnames(G.(zones{1}).nuc), ...
    {'centroids', 'PCA_coeff', 'PCA_latent', 'ellipsoid_center', 'ellipsoid_radii', 'ellipsoid_evecs', 'ellipsoid_v', 'ellipsoid_chi2', 'index', 'PC_range'});

if ~ignore_cells
    for c = 1 : length(characteristics)
        subplot(3, 4, idx);
        idx = idx + 1;
        vals = [];
        for z = 1 : length(zones)
            tmp = G.(zones{z}).cel.(characteristics{c});
            vals(1:length(tmp),z) = tmp;
        end
        vals(vals == 0) = nan;
        [counts, centers] = hist(vals);
        counts = counts ./ repmat(sum(counts,1), size(counts,1), 1);
        bar(centers, counts);
        legend(zones, 'Location', 'Best');
        
        title([regexprep(characteristics{c}, '_', ' '), ' - cells']);
    end
end

for c = 1 : length(characteristics)
    subplot(3, 4, idx);
    idx = idx + 1;
    vals = [];
    for z = 1 : length(zones)
        tmp = G.(zones{z}).nuc.(characteristics{c});
        vals(1:length(tmp),z) = tmp;
    end
    vals(vals == 0) = nan;
    [counts, centers] = hist(vals);
    
    if size(vals,2) == 1
        counts = counts';
        centers = centers';
    end
    
    counts = counts ./ repmat(sum(counts,1), size(counts,1), 1);
    bar(centers, counts);
    legend(zones, 'Location', 'Best');
    ylim([0, 1]);
    
    title([regexprep(characteristics{c}, '_', ' '), ' - nuclei']);
end

if ~ignore_cells
    characteristics = setdiff(fieldnames(G.(zones{1}).inter), {'centroid_shift_physical_coords', 'centroid_shift_pca_coords', 'n_overlapping'});
    for c = 1 : length(characteristics)
        subplot(3, 4, idx);
        idx = idx + 1;
        vals = [];
        for z = 1 : length(zones)
            tmp = G.(zones{z}).inter.(characteristics{c});
            vals(1:length(tmp),z) = tmp;
        end
        vals(vals == 0) = nan;
        [counts, centers] = hist(vals);
        counts = counts ./ repmat(sum(counts,1), size(counts,1), 1);
        bar(centers, counts);
        
        title([regexprep(characteristics{c}, '_', ' '), ' - inter']);
    end

    subplot(3, 4, idx);
    bar(centers, counts);
    ylim([2,3]);
    axis off;
    legend(zones);
end

print('-RGBImage');

% -------------------------------------------------------------------------------------------------------------------------

% now we plot the same data, but with boxplots:

characteristics = setdiff(fieldnames(G.(zones{1}).nuc), ...
    {'centroids', 'PCA_coeff', 'PCA_latent', 'ellipsoid_center', 'ellipsoid_radii', 'ellipsoid_evecs', 'ellipsoid_v', 'ellipsoid_chi2', 'index', 'PC_range'});

fid = figure;
idx = 1;

if ~ignore_cells
    for c = 1 : length(characteristics)
        subplot(3,length(characteristics),idx);
        idx = idx + 1;
        vals = [];
        for z = 1 : length(zones)
            tmp = G.(zones{z}).cel.(characteristics{c});
            vals = [vals; tmp, z*ones(length(tmp),1)];
        end
        boxplot(vals(:,1), vals(:,2));
        set(gca, 'XTickLabel', zones);
        title([regexprep(characteristics{c}, '_', ' '), ' - cells']);
    end
end

for c = 1 : length(characteristics)
    subplot(3,length(characteristics),idx);
    idx = idx + 1;
    vals = [];
    for z = 1 : length(zones)
        tmp = G.(zones{z}).nuc.(characteristics{c});
        vals = [vals; tmp, z*ones(length(tmp),1)];
    end
    boxplot(vals(:,1), vals(:,2));
    set(gca, 'XTickLabel', zones);
    title([regexprep(characteristics{c}, '_', ' '), ' - nuclei']);
end

if ~ignore_cells
    characteristics = setdiff(fieldnames(G.(zones{1}).inter), {'centroid_shift_physical_coords', 'centroid_shift_pca_coords', 'n_overlapping'});
    for c = 1 : length(characteristics)
        subplot(3,4,idx);
        idx = idx + 1;
        vals = [];
        for z = 1 : length(zones)
            tmp = G.(zones{z}).inter.(characteristics{c});
            vals = [vals; tmp, z*ones(length(tmp),1)];
        end
        boxplot(vals(:,1), vals(:,2));
        set(gca, 'XTickLabel', zones);
        title([regexprep(characteristics{c}, '_', ' '), ' - inter']);
    end
end

print('-RGBImage');

% -------------------------------------------------------------------------------------------------------------------------

% now we plot means and standard deviations:

characteristics = setdiff(fieldnames(G.(zones{1}).nuc), ...
    {'centroids', 'PCA_coeff', 'PCA_latent', 'ellipsoid_center', 'ellipsoid_radii', 'ellipsoid_evecs', 'ellipsoid_v', 'ellipsoid_chi2', 'index', 'PC_range'});

fid = figure;
idx = 1;

if ~ignore_cells
    for c = 1 : length(characteristics)
        subplot(3,length(characteristics),idx);
        idx = idx + 1;
        for z = 1 : length(zones)
            tmp = G.(zones{z}).cel.(characteristics{c});
            M(z,1) = mean(tmp);
            SD(z,1) = std(tmp);
        end
        
        bar(M, 'FaceColor', 'none', 'LineWidth', 1, 'BarWidth', 0.5);
        hold on;
        errorbar(M, SD, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12);
        hold off;
        set(gca, 'XTickLabel', zones);
        title([regexprep(characteristics{c}, '_', ' '), ' - cells']);
    end
end

for c = 1 : length(characteristics)
    subplot(3,length(characteristics),idx);
    idx = idx + 1;
    for z = 1 : length(zones)
        tmp = G.(zones{z}).nuc.(characteristics{c});
        M(z,1) = mean(tmp);
        SD(z,1) = std(tmp);
    end
    bar(M, 'FaceColor', 'none', 'LineWidth', 1, 'BarWidth', 0.5);
    hold on;
    errorbar(M, SD, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12);
    hold off;
    set(gca, 'XTickLabel', zones);
    title([regexprep(characteristics{c}, '_', ' '), ' - nuclei']);
end

if ~ignore_cells
    characteristics = setdiff(fieldnames(G.(zones{1}).inter), {'centroid_shift_physical_coords', 'centroid_shift_pca_coords', 'n_overlapping'});
    for c = 1 : length(characteristics)
        subplot(3,4,idx);
        idx = idx + 1;
        for z = 1 : length(zones)
            tmp = G.(zones{z}).inter.(characteristics{c});
            M(z,1) = mean(tmp);
            SD(z,1) = std(tmp);
        end
        bar(M, 'FaceColor', 'none', 'LineWidth', 1, 'BarWidth', 0.5);
        hold on;
        errorbar(M, SD, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 12);
        hold off;
        set(gca, 'XTickLabel', zones);
        title([regexprep(characteristics{c}, '_', ' '), ' - inter']);
    end
end

print('-RGBImage');
