%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load_HEALY1802_MET.m
%
% Load MET data from Healy and combine into one large file.
%
% This script runs San Nguyen's code to parse the met string for all files
% contained within the data directory.
%
% From Luc Rainville
%%

%Strings specifying which variables we want
flds = {'LA','LO','TT','SA','AT','BP','PA','TW','TI', 'LW', 'SW', 'RH'};
dt = 1/24/60;     %Time interval at which we want the data (1 min).

dir_data = '/Volumes/Data/HLY1802/data/met/'; %Note: stored on external harddrive

savefile_local = 'met.mat';
files = dir([dir_data,'18*.MET']);

if exist(savefile_local,'file')
  update_flag = 1;
  load(savefile_local)
  tmp_time = unique(floor([met.time(end):1/24:now now]));
  day_files = zeros(length(files),1)*NaN;
  for k=1:length(files)
    day_files(k) = datenum(files(k).name(1:6),'yymmdd');
  end
  [~,files_to_do]=intersect(day_files,tmp_time);
  files = files(files_to_do);
  old_met = met;
  clear met
else
  update_flag = 0;
end

flag=0;
met = [];

% hb=waitbar(0,'Loading met files');
for ifile = 1:length(files)
  % waitbar(ifile/length(files),hb)
  clear fname tmp
  fname = files(ifile).name;
  disp(['    loading ' fname]);
  tmp = SN_readShipMET([dir_data fname]);
  
  clear tmp2
  tmp2.time = unique(round(tmp.Time/dt))*dt;
  for mm=1:length(flds)
    if isfield(tmp,flds{mm})
      v = tmp.(flds{mm});
      v(v==-99)=NaN;
      tmp2.(flds{mm}) = bindata_AS(tmp.Time, v, tmp2.time);
    else
      tmp2.(flds{mm}) = met.time*NaN;
    end
  end
  
  if ifile==1
    met = tmp2;
    % deal with the README
    met.README{1} = ['time:  Time in Datenum Format, every ',num2str(dt*24*60,'%0.1f'),' min'];
    nk = 0;
    for k=1:length(tmp.README)
      if ~isempty(strcmp(tmp.README{k}(1:2),flds))
        nk = nk+1;
        met.README{nk} = tmp.README{k};
      end
    end    
  else
    for mm=1:length(flds)
      met.(flds{mm}) = cat(1,met.(flds{mm}), tmp2.(flds{mm}));
    end
    met.time = cat(1,met.time, tmp2.time);
  end  
end


if update_flag
  % disp('Updating')
  index_to_keep = 1:find(old_met.time>=met.time(1),1)-1;

  N = length(met.time);
  for mm=1:length(flds)
    tmp = met.(flds{mm});
    if length(tmp)==N
      tmp0 = old_met.(flds{mm});
      met.(flds{mm})=[tmp0(index_to_keep) ; tmp];
    end
  end
  met.time = [old_met.time(index_to_keep) ; met.time];
end

% Quality control
ii = find(met.LO>-100 | met.LA >88 |  met.LA < 50);
met.LO(ii)=NaN;
met.LA(ii)=NaN;

index = find(diff(met.time)>0);
met.time = met.time(index);
for mm=1:length(flds)  
  met.(flds{mm})=met.(flds{mm})(index);
end

save(savefile_local,'met','-v7.3');
