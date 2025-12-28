function [peakH] = peakMeasure(varargin)    
%% ========================================================================
% 
%                          peakMeasure function
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% This code is to automate the peak height measuring process needed for 
% anlayzing raw GC data produced by a gas chromatograph. The code 
% estimates the peak height from the base level. 
% 
% A user should provide the expected arrival time to Electro Capture 
% Detector (ECD) according to gas type to be detected. If not provided, the
% default arrival time is set as:
%           
       expected_arrival_time = 1.8; % default value [min]
% 
% The base level estimation is based on the MATLAB function 'findchangepts' 
% which finds the time point where the root-mean-squre level of the data 
% change [i.e. diff(y_hz)] fluctuates significantly, which indicates that
% the ECD detects current changes. The signal intensity [hz] right before 
% the fluctation is the base level. 
% 
% The code produces a figure showing the base , peak location & height.
% =========================================================================
%% input variables
for iArg = 1:nargin
    if isstruct(varargin{iArg}); dir_ = varargin{iArg}; end
    if ischar(varargin{iArg}); flName = varargin{iArg}; end
    if isnumeric(varargin{iArg}); expected_arrival_time = varargin{iArg}; end
end

%% load GC data 
fidDataFl = fopen([dir_.data flName],'r');
g = textscan(fidDataFl, '%s','delimiter','\n'); fclose(fidDataFl);
g = g{1};

x_time = nan(length(g),1);
y_hz = nan(length(g),1);


for iL = 1:length(g)
    note = textscan(g{iL}, '%s','delimiter','\t' ); note = note{1};
    
    x_time(iL) = str2num(erase(note{1},char(65279)));
    y_hz(iL) = str2num(erase(note{2},char(65279)));
end

unitDataLength = round(numel(y_hz)/100); % 1/100 of total data [1%]
% if not smoothed, too many local minima due to noise in raw data
y_hz_ = smoothdata(y_hz,'gaussian',unitDataLength); 


%% base level & peak 
[~,ind_expected_arrival_time] = min(abs(x_time-expected_arrival_time));

[ind_change] = findchangepts(diff(y_hz(1:ind_expected_arrival_time)),'Statistic','rms');

ind_base_end = min(ind_change);
ind_base_start = ind_base_end-unitDataLength+1;
baseLevel = mean(y_hz(ind_base_start:ind_base_end));

ind_localMin_All = find(islocalmin(y_hz_)); 

ind_localMin_left = ind_localMin_All(ind_localMin_All<ind_expected_arrival_time);
ind_localMin_right = ind_localMin_All(ind_localMin_All>ind_expected_arrival_time);
if isempty(ind_localMin_right); ind_localMin_right = length(x_time); end

ind_curve = ind_localMin_left(end):ind_localMin_right(1); % indices for the target curve
[~,i]=max(y_hz_(ind_curve));

ind_peak = ind_curve(i); % peak index

peakH = y_hz_(ind_peak)-baseLevel;

%% visualization --------------------------
clf; 

xlSpan = (x_time(ind_peak) - x_time(ind_base_start))*1.2;
xl = [-xlSpan xlSpan] + x_time(ind_peak);
yl = [baseLevel y_hz_(ind_peak)] +[-peakH peakH];

subplot(121)
plot(x_time,y_hz,'-','linewidth',1); hold on;

title(sprintf('%s',flName(1:end-4)),'Interpreter','none')

patch([xl flip(xl)],[yl([1 1]) yl([2 2])] ,'k','FaceAlpha',0.01,'edgealpha',0.4)

xlabel('time [min]')
ylabel('hz')

subplot(122)

plot(x_time,y_hz,'-','linewidth',1); hold on;
plot([x_time(ind_base_start) x_time(ind_peak+unitDataLength)],...
    [baseLevel baseLevel],'r:','linewidth',2);
plot([x_time(ind_peak) x_time(ind_peak)],[baseLevel y_hz_(ind_peak)],'c','linewidth',3)
scatter(x_time(ind_peak),y_hz_(ind_peak),100,'r_','linewidth',3)


txt = sprintf('base: %.4f',baseLevel);
text(x_time(ind_peak),baseLevel-peakH*0.1,txt)

txt = sprintf('peak: %.4f',peakH);
text(x_time(ind_peak),baseLevel+peakH*1.1,txt)


title('base - N_2O peak')
xlim(xl)
ylim(yl)

xlabel('time [min]')
ylabel('hz')


fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 6 3];

flNameFig = sprintf([dir_.fig flName(1:end-4) '.png']);

print(fig,flNameFig,'-dpng')

end