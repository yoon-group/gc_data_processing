clear all; clc; clf
%% ========================================================================
% 
%                   Standard N2O gas samples processing
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% This code analyzes only files with a name starting with 'CAL' to avoid
% analyzing non-standard sample files. Naming rule of the standard files
% should be as following:
% 
%       CAL_(triplicate index #)_(stock percentage).CSV
% 
% For example, the first sample among the triplicates of 30% stock mixture 
% gas sample should be named as:
% 
%       CAL_1_30.CSV
% 
% A user should specify the expected arrival time to Electro Capture 
% Detector (ECD) according to gas type to be detected. For N2O, our GC
% works well with:
%           
        expected_arrival_time = 1.8; % [min] 
% 
% The user should also provide manufactured N2O gas stock concentration: 
% 
        stock = 0.97; %[ppm]
%
% Once run, the saved result file will be used to identify new data
% which have not been analyzed, and do the additional anaysis only for the
% new data. 
% 
% =========================================================================

%% directory setting
dir_.archive = sprintf('/home/public/gcData/');

dir_.results = sprintf('%s/processed/',pwd);
if ~exist(dir_.results,'dir'); mkdir(dir_.results); end

dir_.fig = [dir_.results 'fig/STD/']; 
if ~exist(dir_.fig,'dir'); mkdir(dir_.fig); end

dir_.data = sprintf('%sSTD/',dir_.archive);

%% output file name
flNameOut = sprintf('%sSTD.mat',dir_.results);

%% read standard sample file names 
flNameStruct = dir(dir_.data);
nSample = 0;
for il = 3:numel(flNameStruct)
    if strcmpi(flNameStruct(il).name(1:3),'CAL')
        nSample = nSample + 1;
        flNameList{nSample,1} = flNameStruct(il).name;
    end
end

%% Identify unanalyzed files
if exist(flNameOut,'file')
    load(flNameOut,'T_STD');
    flNameList = setdiff(flNameList,T_STD.flNameList);
end

nSample = length(flNameList);
%% result variables 
cnc = nan(nSample,1);
indTriplicate = nan(nSample,1);
peakH = nan(nSample,1);

%% read standard solutions' concentration
for iSample = 1:nSample
    flName = flNameList{iSample}(1:end-4);
    note = textscan(flName, '%s','delimiter','_' ); note = note{1};

    indTriplicate(iSample) = str2double(sprintf('%s',note{2})); % triplicate index
    cnc(iSample) = str2double(sprintf('%s',note{3}))*stock/100; % percentage [ppm]
end

%% Peak measuring
for iSample = 1:nSample
    fprintf('(%d/%d) STD sample %s --> processing -->',iSample,length(flNameList),flNameList{iSample}(1:end-4))
    peakH(iSample) = peakMeasure(dir_,flNameList{iSample},expected_arrival_time);    
    fprintf('  completed \n')
end

%% save results
T_new = table(flNameList,cnc,indTriplicate,peakH);
if ~exist('T_STD','var')
    T_STD = T_new;
elseif isempty(T_new)
    clc
    fprintf("\n No new standard sample data added. \n Add data into: '%s' \n\n",dir_.data)
else
    T_STD = [T_STD; T_new];
end

p = polyfit(T_STD.peakH,T_STD.cnc,1);
slope = p(1);
offset = p(2);
    
% save(flNameOut,'T_STD','slope','offset')

%% calibration curve figures ==============================================
peakH = T_STD.peakH;
cnc = T_STD.cnc;

%  ------- (hz -> N2O)
p = polyfit(peakH,cnc,1);
slope = p(1);
offset = p(2);
   
SSR = sum((peakH*slope + offset - cnc).^2); % sum of squared error
R2 = 1- SSR / sum((cnc - mean(cnc)).^2); % R^2 value
    
clf; 
% subplot(121)
plot(peakH,slope*peakH+offset,'-'); hold on
scatter(peakH,cnc,'o')

ylim([-0.1 1.1])
xl = xlim; yl = ylim;
text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10,sprintf('y = (%.4e) * x + (%.4e)',slope,offset),'FontSize',14)
text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10*2,sprintf('R^2 = %.4f',R2),'FontSize',14)

ylabel('N_2O [ppm]');
xlabel('Peak Height [hz]')
title('hz -> N_2O')
grid on

%  ------- (N2O -> hz)
% p = polyfit(cnc,peakH,1);
% slope = p(1);
% offset = p(2);
% SSR = sum((cnc*slope + offset - peakH).^2); % sum of squared error
% R2 = 1- SSR / sum((peakH - mean(peakH)).^2); % R^2 value
% 
% subplot(122)
% plot(cnc,slope*cnc+offset,'-'); hold on
% scatter(cnc,peakH,'o')
% 
% xl = xlim; yl = ylim;
% text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10,sprintf('y = %.2f * x + %.2f',slope,offset))
% text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10*2,sprintf('R^2 = %.4f',R2))
% 
% grid on
% 
% xlabel('N2O [ppm]');
% ylabel('peak height [hz]')
% title('N_2O -> hz')

set(gca,'FontSize',14)

fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 6 5];

flNameFig = sprintf('%scalibration.png',dir_.fig);

print(fig,flNameFig,'-dpng')
