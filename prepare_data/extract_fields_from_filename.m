function [type, gen_bg, age, bone_name, preg_num, mouse_num, genotype, body_side, postfix] = extract_fields_from_filename(filename)

filename_regexp = [ ...
    '^(Humerus|RadiusUlna|Radius|Ulna|Femur|FemurTibiaFibula|TibiaFibulaFemur|TibiaFibula|Tibia|Fibula|Mandibula|Rib|DistalTibia|ProximalTibia)' ...
    '_(WT|MDG|GDF5|BMP2|BMP4|BMP2BMP4|Normal|OP|Cancerous|MUT|Src|NA|Gli1creERdTA)_' ...
    '(E14|E15|E16|E17|E18|P0|NB|E20|P1|E21|P2|P3|P4|P6|P8|P10|P12|P14|P16|P18|P20|P24|P28|P32|P40|P52|P42|w4|w5|w6|w7|w8|w12|w40|adult|\d+Y|\d\d\dY|NA)_(.*)\.vff$'];

bonename_regexp = '^p(\d+)m(\d+)_(wt|wthom|wthet|het|mut|dh|wthethet|wtwthet|control|mut1)_(L|R|UA|NA|left|right)(.*)';

filename = regexprep(filename, '\.mat', '.vff');

% checking validity of input
if ~ischar(filename), error('Input is not a string'); end

% allocating output
type = [];
gen_bg = [];
age = [];
bone_name = [];
preg_num = [];
mouse_num = [];
genotype = [];
body_side = [];
postfix = [];

% separating the file name
[~, name, ext] = fileparts(filename);
if length(name) >= 6
    if strcmp(name(1:6), '_ROOT_')
        name = name(7:end);
    end
end
rexp = regexp([name, ext], filename_regexp, 'tokens');

% if the filename does not follow the nomenclature
if isempty(rexp), return; end

% passing the output to the output variables
type = rexp{1}{1};
gen_bg = rexp{1}{2};
age = rexp{1}{3};
bone_name = rexp{1}{4};

if regexpi(age, '^\d\d\dY$')
    age = 'adult';
end

if strcmp(gen_bg, 'Normal') || strcmp(gen_bg, 'OP')
    gen_bg = 'WT';
end

rexp = regexp(bone_name, bonename_regexp, 'tokens');

if ~isempty(rexp)
    preg_num = str2double(rexp{1}{1});
    mouse_num = str2double(rexp{1}{2});
    genotype = rexp{1}{3};
    body_side = rexp{1}{4};
    postfix = ['"', rexp{1}{5}, '"'];
else
    preg_num = [];
    mouse_num = [];
    genotype = '';
    body_side = '';
    postfix = '""';
end
