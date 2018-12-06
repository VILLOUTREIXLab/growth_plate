function write_vtk(filename, vertex, face)

vtk_file = fopen(filename, 'wt');
fprintf(vtk_file, '# vtk DataFile Version 3.0\n');
fprintf(vtk_file, 'vtk output\n');
fprintf(vtk_file, 'ASCII\n');
fprintf(vtk_file, 'DATASET POLYDATA\n');
fprintf(vtk_file, ['POINTS ', num2str(size(vertex,1)), ' float\n']);

segments = [0:10000:size(vertex,1), size(vertex,1)];
if segments(end) == segments(end-1)
    segments(end) = [];
end

n_vertex = cell(size(vertex));
for i = 1 : length(segments)-1
    str = [num2str(i), '/', num2str(length(segments)-1)];
    fprintf(str);
    n_vertex(segments(i)+1:segments(i+1),:) = arrayfun(@(x) sprintf('%g', x), vertex(segments(i)+1:segments(i+1),:), 'UniformOutput', false);
    fprintf(repmat('\b', 1, length(str)));
end

vertex = n_vertex;
vertex = vertex';
vertex = vertex(:);

% IDX = 0;
% STR = cast(zeros(1, ceil((length([vertex{:}]) + length(vertex) + length(vertex)*2/9)), 'uint8'), 'char');

STR = cell(length(1 : 9 : length(vertex)),1);
for i = 1 : 9 : length(vertex)
    idcs = i:min(i+8, length(vertex));
    
    str = char(zeros(1, sum(cellfun(@length, vertex(idcs))) + length(idcs) + 2, 'uint8'));
    idx = 0;
    for j = idcs
        str(idx+1:idx+length(vertex{j})+1) = [vertex{j}, ' '];
        idx = idx + length(vertex{j})+1;
    end
    str(end-1:end) = '\n';
    %     STR(IDX+1:IDX+length(str)) = str;
    %     IDX = IDX + length(str);
    
%     fprintf(vtk_file, str);
    STR{(i-1)/9+1,1} = str;
end
fprintf(vtk_file, [STR{:}]);
fprintf(vtk_file, '\n');


fprintf(vtk_file, ['POLYGONS ', num2str(size(face,1)), ' ', num2str(size(face,1)*4), '\n']);

face = face-1;

% IDX = 0;
% len = floor(log10(face))+1;
% len(isinf(len)) = 1;
% len = sum(len,2) + 7;
% STR = cast(zeros(1, sum(len), 'uint8'), 'char');

STR = cell(size(face,1),1);
for i = 1 : size(face,1)
    STR{i,1} = ['3 ', sprintf('%g', face(i,1)), ' ', sprintf('%g', face(i,2)), ' ', sprintf('%g', face(i,3)), ' \n'];
end
fprintf(vtk_file, [STR{:}]);

% fprintf(vtk_file, ['\nCELL_DATA ', num2str(size(face,1)), '\nPOINT_DATA ', num2str(size(vertex,1)), '\n\n']);
fclose(vtk_file);

