function MET = SN_readShipMET(filename)
% SN_READSHIPMET - reads SIO ships MET data files
%
% SN_READSHIPMET(FILENAME) returns a MET structure of variables described
% for MET data files
% 
% SN_READSHIPMET(DIRNAME) reads all *.MET files in the directory DIRNAME
% 
% SN_READSHIPMET({FILE1, FILE2,...}) reads all files indicated
%
% Created 2012/03/30 San Nguyen
% Updated 2012/06/09 San Nguyen - updated to read MET data for most ships
% not just the ones on the R/V Revelle
% Updated 2015/09/02 San Nguyen - update to combine MET data without losing
% any fields added later on
% Updated 2016/10/19 Emily Shroyer- Corrected Readme Information
% Updated 2016/10/20 San Nguyen - Update to read MET files in the same
% folder to preserve maximum number of fields that might have been added
% after the initial recording, not restricting number of fields to defined
% by the first file.
%

% check if it is a single file or a directory and a set of files
if ischar(filename) % dir or file
    switch exist(filename)
        case 2 % if it is a file
            fid = fopen(filename,'r');
            MET = SN_readShipMET(fid);
            fclose(fid);
        case 7 % if it is a directory
            my_METfiles = dir(fullfile(filename, '*.MET'));
            if isempty(my_METfiles)
                MET = [];
                return
            else
                % prepare to read all files
                MET = cell(size(my_METfiles));
                % read the files in the directory
                for i = 1:length(my_METfiles)
                    disp(['reading ' my_METfiles(i).name]);
                    MET{i} = SN_readShipMET(fullfile(filename,my_METfiles(i).name));
                end
                % combine all files into one MET structure
                MET = SN_combineMET(MET{:});
            end
        otherwise
            error('MATLAB:SN_readShipMET:wrongFILENAME','Invalid file specified.');
    end
elseif iscellstr(filename) % cell of files
    % prepare to read all files
    MET = cell(size(filename));
    % read all files
    for i = 1:length(filename)
        disp(['reading ' filename{i}]);
        MET{i} = SN_readShipMET(filename{i});
    end
    % combine all files into one MET structure
    MET = SN_combineMET(MET{:});
else
    
    if (filename<1)
        error('MATLAB:SN_readShipMET:wrongFID','FID is invalid');
    end
    
    MET = SN_readShipMETFile(filename);
    
    return
end

end

% reading MET files through FID
function MET = SN_readShipMETFile(fid)

% make sure the file marker begins at the start
frewind(fid)
fgetl(fid);
str = fgetl(fid);
METdate = floor(datenum(str(7:end),'dd-mmm-yy  HH:MM:SS'));

fgetl(fid);
str = fgetl(fid);

% get headers
str = strrep(str,'#','');
str = strrep(str,'-','_');
% get the column labels
struct_vars = textscan(str,'%s','EndOfLine',' ');
struct_vars = struct_vars{1};

% because the time have different reading format than other columns we have
% to read in the time as string and all other columns as numbers
data_str = char(ones(1,2*length(struct_vars))); % defining textscan pattern
time_ind = NaN;
for i = 1:length(struct_vars)
    if strcmpi(struct_vars{i},'Time')
        data_str((i-1)*2+(1:2)) = '%s';
        time_ind = i;
    else
        data_str((i-1)*2+(1:2)) = '%f';
    end
end

if isnan(time_ind)
    error('MATLAB:SN_readShipMET','We can not find the time in the data file. Please check');
end

% read the whole file again with the textscan pattern to read in files
frewind(fid)
data = textscan(fid,[data_str '%*[^\n]'],'delimiter',' \t','MultipleDelimsAsOne',true,'CommentStyle','#');
data{time_ind} = datenum(data{time_ind},'HHMMSS')-datenum('000000','HHMMSS')+METdate;

%make sure that the time from one day is not mistaken for the other day
time_inv_ind = find(diff(data{time_ind})<0);
data{time_ind}(time_inv_ind) = data{time_ind}(time_inv_ind) - 1;

MET = struct();

for i = 1:length(struct_vars)
    MET.(struct_vars{i}) = data{i};
end

MET.README = {...
    'Time:  Time in Datenum Format';
    'AT:    Air Temperature (deg C)';
    'BP:    Barometric Pressure (mb)';
    'PR:    Precipitation (mm)';
    'RH:    Relative Humidity (%rh)';
    'RT:    Relative Temperature (deg C)';
    'DP:    Dew Point (deg C)';
    'LD:    Long wave radation dome temperature (deg K)';
    'LB:    Long wave radation body temperature (deg K)';
    'LT:    Long wave radation thermopile voltage (volts)';
    'LW:    Long wave radiation (W/m^2)';
    'SW:    Short wave radiation (W/m^2)';
    'PA:    Surface PAR (uE/sec/m^2)';
    'WS:    Relative Wind Speed (m/s)';
    'WD:    Relative Wind Direction (direction wind is coming from)';
    'TW:    True Wind Speed (m/s)';
    'TI:    True Wind Direction (direction wind is coming from)';
    'TT:    TSG Temperature (deg C)';
    'TC:    TSG Conductivity with slope offset  (mS/cm)';
    'SA:    Salinity (psu)';
    'SD:    Sigma t (kg/m^3)';
    'SV:    Sound Speed (m/s)';
    'TG:    TSG Conductivity without slope offset  (mS/cm)';
    'FI:    Flowmeter (lpm)';
    'PS:    Pressure (psi)';
    'OC:    Oxygen Current (ua)';
    'OT:    Oxygen Temperature (deg C)';
    'OX:    Oxygem (ml/l)';
    'OS:    Oxygen Saturation (ml/l)';
    'FL:    Fluorometer (micro gm/l)';
    'VP:    Vertical reference unit pitch (deg)';
    'VR:    Vertical reference unit roll (deg)';
    'VH:    Vertical reference unit heave (m)';
    'VX:    Ships Trim (deg) ';
    'VY:    Ships List (deg)';
    'GY:    Gyro Heading (deg)';
    'MB:    Water Depth (m Center Beam Depth)';
    'BT:    Knudsen Water depth (m)';
    'LA:    NMEA Latitude';
    'LO:    NMEA Longitude';
    'GT:    GPS time of Day (seconds, 0-86400)';
    'CR:    Ships Course Over Ground (deg)';
    'SP:    Ships Speed Over Ground (knots)';
    'ZD:    GPS Date Time (seconds since 00:00:00 01/01/1970';
    'GA:    GPS Altitude (m above/below mean sea level) ';
    'GS:    GPS Status';
    'SH:    Ashtech heading (deg)';
    'SM:    Ashtech pitch (deg)';
    'SR:    Ashtech roll (deg_';
    'ZO:    Winch wire out (m)';
    'ZS:    Winch wire speed (m/min)';
    'ZT:    Wihch wire tension (lbs)';
    'PZ:    Uncontaiminated seawater pump status (on/off)';
    'IP:    CTD depth (m)';
    'IV:    CTD velocity (m/s)';
    'IA:    CTD Altimeter (m)';  
    };
end

function MET = SN_combineMET(varargin)
% SN_combineMET - combines MET data files in MATLAB format that was
% converted using SN_readShipMET
%
% SN_combineMET(MET1,MET2,MET3,...) returns a MET structure of variables described
% for MET data files
%
% Written 2015/09/02 - San Nguyen snguyen@opg1.ucsd.edu

if nargin < 1 
    MET = [];
    return;
end

if nargin == 1
   MET = varargin{1};
   if length(MET) == 1
       return;
   end
   evalstr = 'MET = SN_combineMET(';
   for i=1:(length(MET)-1)
       evalstr = strcat(evalstr, 'MET(', num2str(i), '),');
   end
   evalstr = strcat(evalstr, 'MET(', num2str(length(MET)), '));');
   eval(evalstr);
   return
end

MET_fields = fieldnames(varargin{1});
for i = 2:nargin
    tmp_fields = fieldnames(varargin{i});
    for j = 1:numel(tmp_fields)
        if ~ismember(tmp_fields{j},MET_fields)
            MET_fields{end+1} = tmp_fields{j};
        end
    end
end

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    evalstr = strcat('MET.', MET_fields{i}, '= [');
    for j=1:(nargin-1)
        if ~isfield(varargin{j},(MET_fields{i}))
            varargin{j}.(MET_fields{i}) = NaN(size(varargin{j}.Time));
        end
        evalstr = strcat(evalstr, 'varargin{', num2str(j), '}.', MET_fields{i}, ';');
    end
    if ~isfield(varargin{nargin},(MET_fields{i}))
        varargin{nargin}.(MET_fields{i}) = NaN(size(varargin{nargin}.Time));
    end
    evalstr = strcat(evalstr, 'varargin{', num2str(nargin), '}.', MET_fields{i}, '];');
    eval(evalstr);
end

if (isfield(varargin{1},'README'))
    MET.README = varargin{1}.README;
end

% sort out time
[~, I] = sort(MET.Time);

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    MET.(MET_fields{i}) = MET.(MET_fields{i})(I);
end

%find unique time
[~, I, ~] = unique(MET.Time);

for i=1:(length(MET_fields))
    if strcmpi(MET_fields{i},'README')
        continue;
    end
    MET.(MET_fields{i}) = MET.(MET_fields{i})(I);
end

end