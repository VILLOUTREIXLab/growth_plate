function vffwrite(I, filename)

fid = fopen(filename, 'w', 'b');

% if isfield(I.header, 'original_dt_is_u16')
%     if I.header.original_dt_is_u16 == 1
%         I.img = int16(single(I.img) + single(intmin('int16')));
%         I.data_type = 'int16';
%         I.header.data_type = 'int16';
%         if isfield(I.header, 'th')
%             I.header.th = I.header.th + double(intmin('int16'));
%         end
%     end
% end

if isfield(I.header, 'origin_in_pix')
    I.header.origin = I.header.origin_in_pix;
    I.header = rmfield(I.header, 'origin_in_pix');
    I.header = rmfield(I.header, 'origin_in_mm');
elseif isfield(I.header, 'origin_in_mm')
    I.header.origin = (I.header.origin_in_mm / I.header.spacing(1) * 10) / 10;
    I.header = rmfield(I.header, 'origin_in_mm');
elseif ~isfield(I.header, 'origin')
    I.header.origin = [0 0 0];
end

I.header = restore_header_format(I.header);
fprintf(fid, '%s', I.header);

if iscell(I.img)
    I.img = cellfun(@(slide) slide',I.img,'UniformOutput',false);
    I.img = cellfun(@(slide) slide(:,end:-1:1),I.img,'UniformOutput',false);
    I.img = I.img(end:-1:1);
    for i = 1 : length(I.img)
        fwrite(fid, I.img{i}, 'int16');
    end
else
    %     I.img = ipermute(I.img, [2,1,3]);
    %     I.img = I.img(:, end:-1:1, end:-1:1);
    fwrite(fid, I.img, class(I.img));
end

fclose(fid);

% -------------------------------------------------------------------------

function new_hdr = restore_header_format(hdr)

name = fieldnames(hdr);
val = struct2cell(hdr);

ind = strcmp('hu2mgml', name);
name = name(~ind);
val = val(~ind);
ind = strcmp('mgml2hu', name);
name = name(~ind);
val = val(~ind);
ind = strcmp('Xradia2GE_unq', name);
name = name(~ind);
val = val(~ind);
ind = strcmp('GE2Xradia_unq', name);
name = name(~ind);
val = val(~ind);

is_mis = strcmp(name, 'miscellaneous');

new_hdr = [];
if any(is_mis)
    for i = 1 : length(val{is_mis})
        new_hdr = [new_hdr, val{is_mis}{i}, char(10)]; %#ok<AGROW>
    end
end

for i = 1 : length(name)
    if is_mis(i), continue; end
    if isnumeric(val{i}) || islogical(val{i})
        if size(val{i},1) == 1
            new_hdr = [new_hdr, name{i}, '=', ...
                regexprep(num2str(val{i}(:)'), '\s*', ' '), ';', char(10)]; %#ok<AGROW>
        else
            new_hdr = [new_hdr, name{i}, '=', mat2str(val{i}), ';', char(10)]; %#ok<AGROW>
        end
    else
        try
            new_hdr = [new_hdr, name{i}, '=', val{i}, ';', char(10)]; %#ok<AGROW>
        catch
            keyboard;
        end
    end
end

new_hdr(end+1:end+2) = [char(12), char(10)];

