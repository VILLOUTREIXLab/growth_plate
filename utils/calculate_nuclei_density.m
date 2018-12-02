function [V, X, Y, Z, point_density] = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points)
%{
Input:
    N - an n*3 matrix with the [x,y,z] coordinates of each of the n points.
    spacing -  1*3 array with the [x,y,z] physical size of the voxels (if N is already in physical units, just use [1,1,1]).
    delta_spacing - a positive real value that specifies the frequency of the sampling of the 3D grid in the output image (the density space / the output argument V).
    exclude_boundary_points - a boolean value that allows the user to choose whether to include or exclude boundary vertices to avoid potential noise in the
    physical boundaries of the sample, default is 'false' (i.e. include boundaries)
Output:
    V - the 3D matrix of estimated densities.
    X, Y, Z - the coordinates of the interpolated space V.
Example:
    N = rand(100, 3);
    spacing = [1 1 1];
    delta_spacing = 0.01;
    exclude_boundary_points = false;
    V = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points);
    figure;
    imagesc(V(:,:,50));
    axis image;
    colormap jet;
%}
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
TRI = delaunay(N(:,1), N(:,2), N(:,3));
tetra_vol = zeros(size(TRI,1),1);
for i = 1 : size(TRI,1)
    tetra_vol(i,1) = abs(det([N(TRI(i,:),:)'; ones(1,4)]))/6;
end
% calculating the density at each vertex (i.e. nucleus center), the units are 1/um^3:
point_density = zeros(size(N,1),1);
for i = 1 : size(N,1)
    point_density(i,1) = 4/sum(tetra_vol(any(TRI == i,2)));
end
% in case the user wants to exclude boundary vertices from the analysis:
original_N = N;
if exclude_boundary_points
   
    % to avoid the noise in the physical boundaries of the sample, we remove all points that belong to a tetrahedron that has an external face:
    TR = triangulation(TRI, N);
    [~, points] = freeBoundary(TR);
    % finding the index of each boundary face vertex in the original vertex list:
    point_indices = nan(size(points,1),1);
    for i = 1 : size(points,1)
        point_indices(i,1) = find(points(i,1) == N(:,1) & points(i,2) == N(:,2) & points(i,3) == N(:,3));
    end
    % removing all points that are in the same tetrahedra as in the variable 'points':
    tetrahedra_to_exclude = any(ismember(TRI, point_indices),2);
    points_to_exclude = unique(TRI(tetrahedra_to_exclude,:));
    N(points_to_exclude,:) = [];
    point_density(points_to_exclude,:) = [];
   
end
% setting the grid for the interpolation:
[X,Y,Z] = ndgrid( ...
    min(original_N(:,1)):spacing(1)*delta_spacing:max(original_N(:,1)), ...
    min(original_N(:,2)):spacing(2)*delta_spacing:max(original_N(:,2)), ...
    min(original_N(:,3)):spacing(3)*delta_spacing:max(original_N(:,3)));
% interpolating the densities:
F = scatteredInterpolant(N(:,1), N(:,2), N(:,3), point_density, 'linear', 'none');
V = F(X,Y,Z);
% preparing the XYZ coordinates to be outputed:
X = squeeze(X(:,1,1));
Y = squeeze(Y(1,:,1));
Z = squeeze(Z(1,1,:));

end
