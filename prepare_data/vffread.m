function varargout = vffread(filename, data_type, scale, read_what)

% uint8 - with rescaled default units (rescaling should be by dividing
% with an integer, so that the histogram would still be smooth) - don't
% know yet how...

max_chars_in_header_key = 1000;

if ~exist('read_what','var'), read_what = []; end
if ~exist('filename','var'), filename = []; end
if ~exist('scale','var'), scale = 1; end

% loading the file
if isempty(filename)
    filename = choose_filename;
    if filename == 0
        varargout(1:nargout) = {[]};
        return;
    end
end
fid = fopen(filename,'r','b');

% loading header
I.header = load_header(fid, max_chars_in_header_key);

% here we change from using the field 'origin' to either mm or pixels:
if isfield(I.header, 'origin')
    I.header.origin_in_mm = I.header.origin .* I.header.spacing;
    I.header.origin_in_pix = I.header.origin;
    I.header = rmfield(I.header, 'origin');
else
    I.header.origin_in_mm = [0 0 0];
    I.header.origin_in_pix = [0 0 0];
end

get_data_type = 0;
if ~exist('data_type','var')
    get_data_type = 1;
elseif strcmp(data_type, '')
    get_data_type = 1;
end
if get_data_type
    if ~isfield(I.header, 'data_type')
        if I.header.bits == 8
            I.header.data_type = 'uint8';
        else
            I.header.data_type = 'int16';
        end
    end
    data_type = I.header.data_type;
end

% updating parameters:
I = update_struct_parameters(I, filename, data_type, scale);

if strcmp(read_what, 'header')
    varargout{1} = I.header;
    fclose(fid);
    return;
end

if strcmp(read_what, 'no_img')
    varargout{1} = I;
    fclose(fid);
    return;
end

% loading volume:
% if strcmp(data_type, 'uint16')
%     I.img = read_volume('int16', I, fid, scale);
%     I.img(I.img < 0) = 0;
%     I.img = cast(I.img, 'uint16');
% else
    I.img = read_volume(data_type, I, fid, scale);
% end

%     I.img = I.img + 32768;
%     I.img = I.img - 32767;
    
% if strcmp(I.header.scanner, 'Xradia') && isa(I.img, 'int16')
%     I.img = uint16(double(I.img) - double(intmax('int16')) + double(intmax('uint16')));
%     I.header.data_type = 'uint16';
%     I.data_type = 'uint16';
% end

if isfield(I.header, 'original_dt_is_u16')
    if I.header.original_dt_is_u16 == 1
        I.img = uint16(single(I.img) - single(intmin('int16')));
        I.data_type = 'uint16';
        I.header.data_type = 'uint16';
        if isfield(I.header, 'th')
            I.header.th = I.header.th - double(intmin('int16'));
        end
    end
end

if strcmp(read_what, 'img')
    varargout{1} = I.img;
    fclose(fid);
    return;
end

% closing the 'basta'
fclose(fid);
if nargout == 1
    varargout{1} = I;
else
    varargout{1} = I.img;
    varargout{2} = I.header;
end


function filename = choose_filename

[filename, pathname] = uigetfile(fullfile(cd, '*.vff'), 'Please select the vff file:');
if filename == 0, return; end
filename = fullfile(pathname, filename);


function new_hdr = extract_header_data(hdr, fid)

tokens = regexpi(hdr(1:end-1), '^([^=])+=(.*);', 'tokens');
is_match = ~cellfun(@isempty, tokens);
matches = [tokens{is_match}]';
matches = [matches{:}]';
name = matches(1:2:end);

val = matches(2:2:end);
nval = cellfun(@str2num, regexprep(val, 'image', ''), 'UniformOutput', false);

is_name_field = strcmp(name, 'name');
is_numeric = ~cellfun(@isempty, nval) & ~is_name_field;
val(is_numeric) = nval(is_numeric);

new_hdr.('miscellaneous') = hdr(~is_match);
for i = 1 : length(val)
    name{i} = regexprep(name{i}, '[/\(\)\-]', '_');
    try
        new_hdr.(name{i}) = val{i};
    catch
        fclose(fid);
        error('Problematic header');
    end
end

deleted_field_idx = strcmp(name, 'image_type');
if any(deleted_field_idx)
    new_hdr = rmfield(new_hdr, name{deleted_field_idx});
end

% there are cases where the 'spacing' field is [1 1 1] and there is an
% additional field (not existing in regular cases) named 'elementsize' that
% does contain the true resolution. Here I standarize the format of the
% header:
if isfield(new_hdr, 'elementsize')
    if isequal(new_hdr.spacing, [1 1 1])
        new_hdr.spacing = new_hdr.spacing * new_hdr.elementsize;
        new_hdr = rmfield(new_hdr, 'elementsize');
    elseif isequal(new_hdr.spacing, [1 1])
        new_hdr.spacing = new_hdr.spacing * new_hdr.elementsize;
        new_hdr = rmfield(new_hdr, 'elementsize');
    else % making sure that spacing does not contain important information
        keyboard;
    end
end


function header = load_header(fid, max_chars_in_header_key)

header = {};
while true
    position = ftell(fid);
    header_string = fread(fid, max_chars_in_header_key, '*int8')';
    if ~any(header_string == 10)
        fclose(fid);
        header_string = char(header_string);
        if strcmp(char(header_string(2:9)), 'TicT@cNu')
            error('Corrupted file');
        else
            save('Suspected file', 'header_string');
            error(['Suspected file ', header_string]);
        end
    else
        fseek(fid, position, -1);
    end
    
    header{end+1,1} = fgetl(fid); %#ok<AGROW>
    
    if ~regexp(header{end}, '^[\"\%\-\.\/0123456789\:\;\=ABCDEFGHIJKLMNOPQRSTUVWXYZ\\\_abcdefghijklmnopqrstuvwxyz]+$')
        keyboard;
    end
    
    if header{end} == -1
        fclose(fid);
        error('Corrupted image file');
    end
    
    if strcmp(header{end}, char(12))
        break;
    end
end
header = extract_header_data(header, fid);

% % completing missing fields based on the xls preferences file:
% fn = fieldnames(P.vff);
% for i = 1 : length(fn)
%     if ~isfield(header, fn{i})
%         header.(fn{i}) = P.vff.(fn{i});
%     end
% end

% if isfield(header, 'boneHU')
%     header.hu2mgml = @(x) single(x) * P.calibration_val  / (header.boneHU - header.water);
%     header.mgml2hu = @(x) single(x) * (header.boneHU - header.water) / P.calibration_val;
% else
%     header.hu2mgml = @(x) x;
%     header.mgml2hu = @(x) x;
% end

header = orderfields(header);
if length(header.size) == 2 % it's a MIP
    header.size(3) = 1;
end


function I = update_struct_parameters(I, filename, data_type, scale)

I.filename = filename;
I.size = floor(I.header.size / scale);
I.header.data_type = data_type;
I.data_type = data_type;
I.scale = scale;
I.spacing = I.header.spacing * scale;
I.origin = [repmat(I.scale-1,1,3); mod(I.header.size, I.scale)];

% extracting the id of the scan in which the bone was scanned:
I.scan_id = '';
if isfield(I.header, 'cmdLine')
    tok = regexp(I.header.cmdLine, ' -i (\S+) ', 'tokens');
    if ~isempty(tok)
        I.scan_id = tok{1}{1};
        if length(tok) > 1 || length(tok{1}) > 1
            keyboard;
        end
    end
end

[I.type, I.gen_bg, I.age, I.bone_name, I.preg_num, I.mouse_num, I.genotype, I.body_side] = extract_fields_from_filename(filename);

% if isempty(I.type)
%     I.SE_shape = [];
%     I.SE_size = [];
%     I.th_sd = [];
% else
%     I.SE_shape = set_SE_shape(I.age);
%     I.SE_size = set_SE_size(I.age);
%     I.th_sd = set_th_sd(I.type, I.age);
% end

% in case the scanner of the image is not given within the header:
if ~isfield(I.header, 'scanner')
    I.header.scanner = '';
end


function img = read_single(I, fid, scale)

img = zeros(I.size(1),I.size(2),I.size(3),'single');
seg = zeros(I.header.size(1), I.header.size(2), 'single');
seg_size = prod(I.header.size(1:2));

if scale == 1
    img = single(fread(fid, inf, 'int16'));
    img = I.header.hu2mgml(img);
else
    for i = 1 : I.size(3)
        fseek(fid, seg_size*2*(scale-1), 'cof');
        seg(:) = single(fread(fid, seg_size, 'int16'));
        seg = I.header.hu2mgml(seg);
        img(:,:,i) = seg(scale:scale:end,scale:scale:end);
    end
end


function img = read_int16(I, fid, scale)

img = zeros(I.size(1),I.size(2),I.size(3),'int16');
seg = zeros(I.header.size(1), I.header.size(2), 'int16');
seg_size = prod(I.header.size(1:2));

if scale == 1
    fseek(fid, seg_size*2*(scale-1), 'cof');
    img(:) = cast(fread(fid, seg_size*I.size(3), 'int16'), 'int16');
else
    for i = 1 : I.size(3)
        fseek(fid, seg_size*2*(scale-1), 'cof');
        tmp = cast(fread(fid, seg_size, 'int16'), 'int16');
        if numel(tmp) ~= numel(seg)
            keyboard;
        end
        seg(:) = tmp;
        img(:,:,i) = seg(scale:scale:end,scale:scale:end);
    end
end


function img = read_uint16(I, fid, scale)

img = zeros(I.size(1),I.size(2),I.size(3),'uint16');

seg = zeros(I.header.size(1), I.header.size(2), 'uint16');
seg_size = prod(I.header.size(1:2));
for i = 1 : I.size(3)
    fseek(fid, seg_size*2*(scale-1), 'cof');
    try
        seg(:) = cast(fread(fid, seg_size, 'uint16'), 'uint16');
    catch
        keyboard;
    end
    img(:,:,i) = seg(scale:scale:end,scale:scale:end);
end


function img = read_uint8(I, fid, scale)

img = zeros(I.size(1),I.size(2),I.size(3),'uint8');
seg = zeros(I.header.size(1), I.header.size(2), 'single');
seg_size = prod(I.header.size(1:2));

for i = 1 : I.size(3)
    fseek(fid, seg_size*(scale-1), 'cof');
    try
        if i == I.size(3)
            tmp = fread(fid, inf, '*uint8');
            if length(tmp) ~= seg_size
                disp('Last slide misses bytes');
                seg(1:length(tmp)) = tmp;
                seg(length(tmp)+1:end) = 0;
            end
        else
            seg(:) = fread(fid, seg_size, '*uint8');
        end
    catch
        keyboard;
    end
    img(:,:,i) = seg(scale:scale:end,scale:scale:end);
end

% for i = 1 : I.size(3)
%     fseek(fid, seg_size*2*(scale-1), 'cof');
%     seg(:) = cast(fread(fid, seg_size, 'uint8'), 'single');
%     tmp = seg(scale:scale:end,scale:scale:end);
%     tmp = I.header.hu2mgml(tmp) * P.uint8_coeff;
%     tmp(tmp < intmin('uint8')) = intmin('uint8');
%     tmp(tmp > intmax('uint8')) = intmax('uint8');
%     img(:,:,i) = cast(tmp, 'uint8');
% end


function img = read_volume(data_type, I, fid, scale)

% reading volume:
% (god knows why, sometimes it misses the first byte...)

fseek(fid, -2*prod(I.header.size), 'eof');
if strcmp(data_type,'single');
    img = read_single(I, fid, scale);
elseif strcmp(data_type,'int16')
    img = read_int16(I, fid, scale);
elseif strcmp(data_type,'uint8')
    img = read_uint8(I, fid, scale);
elseif strcmp(data_type,'uint16')
    img = read_uint16(I, fid, scale);
end

img = reshape(img, I.size);

% The following paragraph is a brief background on vff file format copied
% from: "http://www.cns.nyu.edu/~jonas/doc/SurfRelax-HOWTO-5.html", hope
% it's correct.
%
% "The VFF format is a modified version of the (now defunct) SunVision volume
% file format (vff). VFF files have a plain text header that describes image
% dimensions, voxel sizes, and origin (the latter coded in a rather non-intuitve
% manner), offset and scale factor. The header is arbitrary length and can
% contain any combination of keyword/value pairs. The header is terminated by
% a CTRL-L character (page break), and is followed directly by binary image
% data stored as unsigned shorts or unsigned chars. The voxel values are stored
% in the following order: y changes fastest, then z, then x. The first read
% voxel is in the upper lefthand anterior corner of the image volume. VFF files
% are primarily used for storing 1D data, but some of the intermediate files
% used by SurfRelax are stored as 3D VFF files."
