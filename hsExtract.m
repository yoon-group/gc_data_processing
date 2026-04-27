clear all; clc; clf
%% ========================================================================
% 
%                   Head-space gas GC data processing
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% This code estimates peak heights of head-space sample GC scans.
% 
% Naming rule of the sample data files should be as following:
% 
%  (text 1)_(text 2)_(text 3)_(text 4)_(text 5).CSV
% 
%  where: 
%       text 1: location (e.g., RR, W3, W4)
%       text 2: sampling volume (e.g., 10ML, 25ML, 72ML, 158ML)
%       text 3: inert gas (e.g., N2, HE)
%       text 4: shaking time (e.g. S3, S5, S8).
%       text 3: tripcliate index (e.g. 03, 05, 08) 
%       text 5: temperature (e.g. T24.9)
% 
% example:
%       RR_72ML_HE_S3_01_T22.9.CSV
% 
% A user should specify the expected arrival time to Electro Capture 
% Detector (ECD) according to gas type to be detected. For N2O, our GC
% works with:
%           
    expcted_arrival_time_min  = 1.89;% [min] 
    expcted_arrival_time = 1.9; % [min] 
%
% Once run, the saved result file will be used to identify new data
% which have not been analyzed, and do the additional anaysis only for the
% new data. 
% =========================================================================
%% directory setting
dir_.data = sprintf('/home/public/gcData/hsExtract/Apr2026/'); % data folder

dir_.results = sprintf('%s/processed_hsExtract/Apr2026/',pwd); % where to save results
if ~exist(dir_.results,'dir'); mkdir(dir_.results); end

dir_.fig = [dir_.results 'fig/']; % where to save figures
if ~exist(dir_.fig,'dir'); mkdir(dir_.fig); end



%% output file name
flNameOut = sprintf('%ssample.mat',dir_.results);

%% read sample file names
flNameStruct = dir(dir_.data);
nSample = 0;
for il = 3:numel(flNameStruct)
    if strcmpi(flNameStruct(il).name(end-2:end),'CSV')
        nSample = nSample + 1;
        flNameList{nSample,1} = flNameStruct(il).name;
    end
end

%% Identify unanalyzed files
% if exist(flNameOut,'file')
%     load(flNameOut,'T_sample');
%     flNameList = setdiff(flNameList,T_sample.flNameList);
% end
nSample = length(flNameList);

%% peak measuring
% expcted_arrival_time_min  = 1.875;% [min] 
peakH = nan(nSample,1);
peakT = nan(nSample,1);
for iSample = 1:length(flNameList)
    fprintf('(%d/%d) %s --> processing -->',iSample,length(flNameList),flNameList{iSample}(1:end-4))
    [peakH(iSample), peakT(iSample)] =peakMeasure(flNameList{iSample},...
                                        expcted_arrival_time_min,...
                                        expcted_arrival_time,dir_);    
    fprintf('  completed \n')
end

%% add newly processed data & save
if nSample
    % T_new = table(flNameList,inertGas,location,bottleSize,shakeTime,ind_triplicate,...
    %     tempK,peakH);
    T_new = table(flNameList,peakT,peakH);
    if ~exist('T_sample','var')
        T_sample = T_new;    
    else
        T_sample = [T_sample; T_new];
    end
    save(flNameOut,'T_sample')
else
    fprintf("\n No new gas sample data added. \n Add new data into: '%s'\n",...
            dir_.data)
end

%%
tmp = T_sample.peakT;
min(tmp(tmp>1.85))