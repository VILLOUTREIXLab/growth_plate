function [files, dirname] = get_file_name(file_type, message, MultiSelect)

% loading the name of the file that has the last loaded directory:
last_dir_filename = get_last_dir_filename;

% loading the path of the dir used in the last operation of the script (unless it's the first use):
if exist(last_dir_filename, 'file')
    load(last_dir_filename);
else
    last_dir = 'C:\Users\colla\Google Drive\Collar project DB\Preprocessed';
end

% asking the user to select files:
if exist('MultiSelect', 'var')
    if isnumeric(MultiSelect)
        if MultiSelect ~= 0
            MultiSelect = 'on';
        else
            MultiSelect = 'off';
        end
    end
else
    MultiSelect = 'on';
end

% setting the default file type for all types:
if ~exist('file_type', 'var')
    file_type = '*.*';
end

% setting the default message:
if ~exist('message', 'var')
    message = 'Please select the files:';
end

% asking the user to select the files:
[files, dirname] = uigetfile(file_type, message, last_dir, 'MultiSelect', MultiSelect);

% validating the response of the user:
if isempty(files) || isnumeric(files)
    error('No files were selected, aborting script.');
end

% in case a single file has been selected, we put it inside a cell:
if isa(files, 'char')
    files = {files};
end

files = files(:);

% saving the path of the dir for the next call to this function:
last_dir = dirname; %#ok<NASGU>
save(last_dir_filename, 'last_dir');
