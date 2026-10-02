clear all; clc; clf
%% ========================================================================
% 
%                   Head-space gas sample processing
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% This code estimates N2O concentration of a gas sample collected from
% head-space (HS) gas extraction method. 
% 
% First run hsExtract.m to collect peak height values from all GC scans.
% 
% A user should provide the value of volumes of aqueous phase (Vaq) and gas 
% phase (Vg) in the HS protocol. If not provided this code assumes as:
% 
    % Vaq : Vg = 9 : 1 
% 
% This code uses the following empirical constants for calculating the 
% Henry's law for N2O (Sander et a., 2006, JPL-Publ-06-2, NASA Technical 
% Reports Server, https://ntrs.nasa.gov/citations/20090033862)
% 
    A = -148.1; B = 8610; C = 20.266; 
% 
% These constants describe solubility equilibrium between the gas phase and 
% the aqueous phase as a function of temperature. The specific formula used
% is:
% 
%       logK0 = A + B / tempK + C * log(tempK)
%
% where tempK is the temperature in Kelvin. The values of A, B, and C are
% specific to N2O.  
% 
% The universal gas constant is given as:
% 
    R = 0.0821; %[L*atm/K/mol]
% 
% Atomic mass units of nitrogen N, O, and N2O [amu] are given as:
% 
    amuN = 14.0067; 
    amuO = 16;
    amuN2O = 2*amuN + amuO;
% =========================================================================
%% load standard sample data and calibration curves
flNameSTD = sprintf('%s/processed_STD/STD.mat',pwd);
load(flNameSTD,'T_STD','slope','offset')

%% load sample data
flNameSample = sprintf('%s/processed_hsExtract/Apr2026/sample.mat',pwd);
load(flNameSample,'T_sample')

T_sample_1 = T_sample; clear T_sample


flNameSample = sprintf('%s/processed_hsExtract/June2026/sample.mat',pwd);
load(flNameSample,'T_sample')
T_sample_2 = T_sample; clear T_sample

T_sample = [T_sample_1; T_sample_2(1:end-2,:)];


%%
flNameList = T_sample.flNameList;
nSample = numel(flNameList);
for iSample = 1:nSample

    note = textscan(flNameList{iSample}(1:end-4), '%s','delimiter','_' ); 
    note = note{1};

    order(iSample,1) = str2num(note{1});
    location{iSample,1} = note{2};
    bottleSize(iSample,1) = str2num(note{3});
    ratio{iSample,1} = note{4};
    shakingTime{iSample,1} = note{5};
    standingTime{iSample,1} = note{6};
    triplicate{iSample,1} = note{7};
    tempK(iSample,1) = str2num(erase(note{8},'T')) + 273.15;
    

    
end
%%
T_sample = addvars(T_sample,order,location,bottleSize,ratio,...
                            shakingTime,standingTime,tempK,triplicate);
[~,ind]=sort(order);

T_sample = T_sample(ind,:);
%% weird sample 
% % #97 from checking fitting figure;
% 97 220 238 268

idx = find(ismember(T_sample.order, [ 97 220 238  ...
    293 294 ...
    295 296 ...
    290 292 ]));
T_sample(idx,:) = [];


flNameList = T_sample.flNameList;
%%
location = {'RR','W3','W4'};
bottleSize = {'158', '72', '25', '10'};
ratio = {'19','441'};
shakingTime = {'S5', 'S13'};
standingTime = {'00', '63'};
triplicate = {'01','02','03'};

iSample = 0;
for iRat = 1:numel(ratio)   
for iLoc = 1:numel(location)
    for iSha= 1:numel(shakingTime)

    
    for iBot = 1:numel(bottleSize)
    
    for iTri = 1:numel(triplicate)
    for iSta = 1:numel(standingTime)
        iSample = iSample + 1;
        flNameList_all{iSample,1}  = sprintf('%d_%s_%s_%s_%s_%s_%s_T', iSample,...
            location{iLoc},...
            bottleSize{iBot},...
            ratio{iRat},...
            shakingTime{iSha},...
            standingTime{iSta},...
            triplicate{iTri});

        flNameList_all{iSample,2} = sum((contains(flNameList,flNameList_all{iSample,1})));
    
    end    
    end
    end
    end
end
end 
%%
fprintf('========= lost samples:\n')
indMissing = find(vertcat(flNameList_all{:,2})==0);
for ind = 1:numel(indMissing)
    fprintf('\t%s\n',flNameList_all{indMissing(ind)})
end


%% ======================================================================== 
tempK = T_sample.tempK;
peakH = T_sample.peakH;

bottleVol = T_sample.bottleSize/1000; %[L]

isRatio_19 = strcmp(T_sample.ratio,'19');
isRatio_441 = strcmp(T_sample.ratio,'441');

nSample = size(T_sample,2);
Vaq = nan(nSample,1); Vg = nan(nSample,1);

Vaq(isRatio_19) = bottleVol(isRatio_19)*9/10; %[L] 
Vg(isRatio_19) =  bottleVol(isRatio_19)*1/10; %[L] 

Vaq(isRatio_441) = bottleVol(isRatio_441)*1/5.4; %[L] 
Vg(isRatio_441) =  bottleVol(isRatio_441)*4.4/5.4; %[L] 


valveVol = 0.3/1000; %[L]
% n2o_ppm_valve = 0/1000; %[ppm]

logK0 = A + (B./tempK) + C*(log(tempK));
K0 = exp(logK0); % Henry's law constant [mol/L/atm]
    
    % Not valve corrected -----------------------------------------
    n2o_ppm_uncorrected = (peakH*slope + offset);
    
    % Valve corrected -----------------------------------------
    n2o_ppm = (n2o_ppm_uncorrected .* (Vg+valveVol) )./(Vg);
    
    n2o_atmUnit = n2o_ppm * 1e-6; % Fractional concentration of N2O in the gas 
                              % phase. By mulitiplying the total pressure 
                              % of 1 atmospheres (atm), the value also 
                              % expresses a partial pressure [atm].
                              % (Henry's law: naq/Vaq = KO * n2o_atm)

    
    ng = (n2o_atmUnit .* Vg)./(R*tempK); % N2O moles in gas phase (from ideal gas law)
    naq = K0.*n2o_atmUnit.*Vaq; % dissolved N2O moles in aqueous phase (from Henry's law)
    
    n2o_mol_L = (ng+naq)./Vaq; % Total N2O concentration in mol/L
    n2o_molN_L = n2o_mol_L * 2;
    n2o_gN_L = n2o_molN_L * amuN;
    n2o_mgN_L = n2o_gN_L * 1e+3;
    n2o_ugN_L = n2o_mgN_L * 1e+3;

    % atmospheric N2O -----------------------------------------
    % logK0_atm = A + (B./mean(tempK)) + C*(log(mean(tempK)));
    % K0_atm = exp(logK0_atm); % Henry's law constant [mol/L/atm]
    % 
    % n2o_atm_atmUnit = n2o_ppm_valve * 1e-6;
    % ng_atm = (n2o_atm_atmUnit .* Vg)./(R*mean(tempK)); % N2O moles in gas phase (from ideal gas law)
    % naq_atm = K0_atm.*n2o_atm_atmUnit.*Vaq; % dissolved N2O moles in aqueous phase (from Henry's law)
    % 
    % n2o_atm_mol_L = (ng_atm+naq_atm)./Vaq; % Total N2O concentration in mol/L
    % n2o_atm_molN_L = n2o_atm_mol_L * 2;
    % n2o_atm_gN_L = n2o_atm_molN_L * amuN;
    % n2o_atm_mgN_L = n2o_atm_gN_L * 1e+3;
    % n2o_atm_ugN_L = n2o_atm_mgN_L * 1e+3;

% into data table
% T_sample = removevars(T_sample,'n2o_ppm_measured');
% T_sample.n2o_ppm_uncorrected = n2o_ppm_uncorrected;
T_sample.n2o_ppm = n2o_ppm;
T_sample.n2o_ugN_L = n2o_ugN_L;
% T_sample.n2o_atm_ugN_L = n2o_atm_ugN_L;



%% data arrangements



n2o_avg_loc_vol_ratio_sht_stt = nan(numel(location),numel(bottleSize),...
                      numel(ratio),numel(shakingTime),numel(standingTime));
n2o_std_loc_vol_ratio_sht_stt = n2o_avg_loc_vol_ratio_sht_stt;
sample_count = n2o_avg_loc_vol_ratio_sht_stt;

% for iLoc = 3%1:numel(location)
% for iVol = 2%1:numel(bottleSize)
% for iRat = 2%1:numel(ratio)
% for iSht = 2%1:numel(shakingTime)
% for iStt = 2%1:numel(standingTime)
for iLoc = 1:numel(location)
for iVol = 1:numel(bottleSize)
for iRat = 1:numel(ratio)
for iSht = 1:numel(shakingTime)
for iStt = 1:numel(standingTime)
    tmp = T_sample(strcmp(T_sample.location,location{iLoc})... 
           & T_sample.bottleSize==str2num(bottleSize{iVol}) ...
           & strcmp(T_sample.ratio,ratio{iRat}) ...
           & strcmp(T_sample.shakingTime,shakingTime{iSht}) ...
           & strcmp(T_sample.standingTime,standingTime{iStt}) ...
           ,:).n2o_ugN_L;
   
    % T_sample(strcmp(T_sample.location,location{iLoc})... 
    %        & T_sample.bottleSize==str2num(bottleSize{iVol}) ...
    %        & strcmp(T_sample.ratio,ratio{iRat}) ...
    %        & strcmp(T_sample.shakingTime,shakingTime{iSht}) ...
    %                   & strcmp(T_sample.standingTime,standingTime{iStt}) ...
    %        ,:).order
    % 

    n2o_avg_loc_vol_ratio_sht_stt(iLoc,iVol,iRat,iSht,iStt) = mean(tmp);
    n2o_std_loc_vol_ratio_sht_stt(iLoc,iVol,iRat,iSht,iStt) = std(tmp);
    sample_count(iLoc,iVol,iRat,iSht,iStt) = numel(tmp);
    % n2o_atm_avg_vol_time(iGas,iLoc,iVol,iTime) = mean(n2o_atm_ugN_L( strcmp(T_sample.inertGas,inertGas{iGas})... 
    %        & strcmp(T_sample.location,location{iLoc}) ...
    %        & T_sample.bottleSize==bottleSize(iVol) ...    %        & T_sample.shakeTime==shakeTime(iTime) ));
    
    
end
end
end
end
end
% writetable(T_sample,'HS_extract.xlsx')


%% paper Figure 3
% (location),(bottleSize),(ratio),(shakingTime),(standingTime)
fs = 20;
clf; clc;

location_txt{1} = 'River Rise';
location_txt{2} = 'Well 3';
location_txt{3} = 'Well 4';


% =======================================
% iCase = 0;
% for ratTxt = {'4.4:1','1:9'}
%     for shtTxt = {'5','13'}
%         for stdTxt = {'0','63'}
%             iCase = iCase + 1;
%             lgndTxt{iCase} = sprintf('%7s %4s (shk) %4s (stnd)',ratTxt{1},shtTxt{1},stdTxt{1});            
%         end
%     end
% end

clear lgndTxt
ratTxt = {'4.4:1','1:9'};
shtTxt = {'13','5'};
stdTxt = {'63','0'};
lgndTxt{1} = sprintf('%7s (ratio) %4s (shaking) %4s (standing)',ratTxt{1},shtTxt{1},stdTxt{1});            
% lgndTxt{2} = sprintf('%7s %4s (shk) %4s (stnd)',ratTxt{1},shtTxt{1},stdTxt{2});            
lgndTxt{2} = sprintf('%7s (ratio) %4s (shaking) %4s (standing)',ratTxt{2},shtTxt{2},stdTxt{2});            
lgndTxt{3} = sprintf('%7s (ratio) %4s (shaking) %4s (standing)',ratTxt{2},shtTxt{2},stdTxt{1});            
% =======================================

plotMat = reshape(1:90,18,5)';
locMat{1} = plotMat(:,1:6);
locMat{2} = plotMat(:,7:12);
locMat{3} = plotMat(:,13:18);


for iLoc = 1:3
    
data_19_Sh5_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,1,1));
data_19_Sh5_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,1,1));
sample_count_19_Sh5_St0= squeeze(sample_count(iLoc,:,1,1,1));

data_19_Sh5_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,1,2));
data_19_Sh5_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,1,2));
sample_count_19_Sh5_St63= squeeze(sample_count(iLoc,:,1,1,2));

data_19_Sh13_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,2,1));
data_19_Sh13_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,2,1));

data_19_Sh13_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,2,2));
data_19_Sh13_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,2,2));



data_441_Sh5_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,1,1));
data_441_Sh5_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,1,1));

data_441_Sh5_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,1,2));
data_441_Sh5_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,1,2));

data_441_Sh13_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,2,1));
data_441_Sh13_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,2,1));

data_441_Sh13_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,2,2));
data_441_Sh13_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,2,2));
sample_count_441_Sh13_St63= squeeze(sample_count(iLoc,:,2,2,2));



% clf
subplot(6,18,locMat{iLoc}(:))
% errorbar(flip([10:10:40])-2,data_441_Sh5_St0,data_441_Sh5_St0_std,'ro:','markersize',10,'linewidth',2); hold on
% errorbar(flip([10:10:40])-1,data_441_Sh5_St63,data_441_Sh5_St63_std,'go:','markersize',10,'linewidth',2); 
% errorbar(flip([10:10:40])+1,data_441_Sh13_St0,data_441_Sh13_St0_std,'cv-.','markersize',10,'linewidth',2); 
errorbar(flip([10:10:40]),data_441_Sh13_St63, data_441_Sh13_St63_std,'ks-','markersize',10,'linewidth',2); hold on


errorbar(flip([10:10:40]),data_19_Sh5_St0,data_19_Sh5_St0_std,'ro:','markersize',10,'linewidth',2); hold on
% errorbar(flip([10:10:40]),data_19_Sh13_St63,data_19_Sh13_St63_std,'g^--','markersize',10,'linewidth',2); 
errorbar(flip([10:10:40]),data_19_Sh5_St63,data_19_Sh5_St63_std,'b^--','markersize',10,'linewidth',2); 
% % errorbar([10:10:40]+0,data_19_Sh13_St0,data_19_Sh13_St0_std,'ks-','linewidth',2); 


xPos = flip([10:10:40]);
% for i = 1:length(xPos)
%     % % (4.4:1, 13, 63)
%     text(xPos(i)+1, ...
%          data_441_Sh13_St63(i) + 0.2, ...
%          num2str(sample_count_441_Sh13_St63(i)), ...
%          'Color','k', ...
%          'FontSize',fs, ...
%          'HorizontalAlignment','center');
% 
%     if iLoc == 1
%         text(25, mean(data_441_Sh13_St63(2:end))+0.4, ...
%              num2str(mean(data_441_Sh13_St63(2:end))), ...
%              'Color','k', ...
%              'FontSize',fs, ...
%              'HorizontalAlignment','center');
%     else
%         text(25, mean(data_441_Sh13_St63(:))+0.4, ...
%              num2str(mean(data_441_Sh13_St63(:))), ...
%              'Color','k', ...
%              'FontSize',fs, ...
%              'HorizontalAlignment','center');
%     end
%     % % (1:9, 5, 0)
%     text(xPos(i)-1, ...
%          data_19_Sh5_St0(i) - 0.2, ...
%          num2str(sample_count_19_Sh5_St0(i)), ...
%          'Color','r', ...
%          'FontSize',fs, ...
%          'HorizontalAlignment','center');
% 
%     text(25, mean(data_19_Sh5_St0(:))-0.4, ...
%          num2str(mean(data_19_Sh5_St0(:))), ...
%          'Color','r', ...
%          'FontSize',fs, ...
%          'HorizontalAlignment','center');
% 
%     % % (1:9, 5, 63)
%     text(xPos(i)+1, ...
%          data_19_Sh5_St63(i) - 0.2, ...
%          num2str(sample_count_19_Sh5_St63(i)), ...
%          'Color','b', ...
%          'FontSize',fs, ...
%          'HorizontalAlignment','center');
% 
%     text(25, mean(data_19_Sh5_St63(:))-0.6, ...
%          num2str(mean(data_19_Sh5_St63(:))), ...
%          'Color','b', ...
%          'FontSize',fs, ...
%          'HorizontalAlignment','center');
% 
% end




set(gca, 'XTick', [10 20 30 40], ...
         'XTickLabel', [],...
         'fontSize',fs)
xlim([5 45])
% ylim([0 4])
% xlabel('bottle size [mL]')
if iLoc ==1; ylabel('N_2O [\mug N/L]'); end
title([location_txt{iLoc}])


if iLoc ~= 1
    set(gca,'YTickLabel',[])
end

% if iLoc == 3
%     legend(lgndTxt{1},lgndTxt{2},lgndTxt{3},lgndTxt{4},'interpreter','none','location','north')
%     legend box off
% end
if iLoc ~= 1
    set(gca,'YTickLabel',[])
end




% subplot(6,18,locMat{iLoc,2}(:))

set(gca, 'XTick', [10 20 30 40], ...
         'XTickLabel', flip(bottleSize),...
         'fontSize',fs)
xlim([5 45])
ylim([-0.1 4])

xlabel('bottle size [mL]')
if iLoc ==1; ylabel('N_2O [\mug N/L]'); end
% title([location{iLoc}])


if iLoc == 2
    legend(lgndTxt,'location','north','fontsize',fs)
    legend box off
end
% 
% if iLoc == 3
%     legend(lgndTxt{5},lgndTxt{6},lgndTxt{7},lgndTxt{8},'interpreter','none','location','north')
%     legend box off
% end
if iLoc ~= 1
    set(gca,'YTickLabel',[])
end


% set(gca,'XTicks', 1:4)

end


fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 24 8];

flNameFig = sprintf('%s/processed_hsExtract/Apr2026/fig_bottleSize_err.png',pwd);


print(fig,flNameFig,'-dpng')

% boxplot([data_19_Sh5_St0(:) data_441_Sh5_St0(:)])% dataW3_19(:) dataW3_441(:) dataW4_441(:) dataW4_19(:)])
% hold on;

%% error bar 
% % (location),(bottleSize),(ratio),(shakingTime),(standingTime)
% 
% clf; clc;
% 
% iCase = 0;
% for ratTxt = {'1:9', '4.4:1'}
%     for shtTxt = {'5','13'}
%         for stdTxt = {'0','63'}
%             iCase = iCase + 1;
%             lgndTxt{iCase} = sprintf('%7s %4s (shk) %4s (stnd)',ratTxt{1},shtTxt{1},stdTxt{1});
%         end
%     end
% end
% 
% plotMat = reshape(1:108,18,6)';
% locMat{1,1} = plotMat(1:3,1:6);
% locMat{1,2} = plotMat(4:6,1:6);
% locMat{2,1} = plotMat(1:3,7:12);
% locMat{2,2} = plotMat(4:6,7:12);
% locMat{3,1} = plotMat(1:3,13:18);
% locMat{3,2} = plotMat(4:6,13:18);
% 
% for iLoc = 1:3
% 
% data_19_Sh5_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,1,1));
% data_19_Sh5_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,1,1));
% 
% data_19_Sh5_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,1,2));
% data_19_Sh5_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,1,2));
% 
% data_19_Sh13_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,2,1));
% data_19_Sh13_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,2,1));
% 
% data_19_Sh13_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,1,2,2));
% data_19_Sh13_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,1,2,2));
% 
% 
% 
% data_441_Sh5_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,1,1));
% data_441_Sh5_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,1,1));
% 
% data_441_Sh5_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,1,2));
% data_441_Sh5_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,1,2));
% 
% data_441_Sh13_St0= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,2,1));
% data_441_Sh13_St0_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,2,1));
% 
% data_441_Sh13_St63= squeeze(n2o_avg_loc_vol_ratio_sht_stt(iLoc,:,2,2,2));
% data_441_Sh13_St63_std= squeeze(n2o_std_loc_vol_ratio_sht_stt(iLoc,:,2,2,2));
% 
% % clf
% subplot(6,18,locMat{iLoc,1}(:))
% errorbar([10:10:40]-2,data_19_Sh5_St0,data_19_Sh5_St0_std,'rs-','linewidth',2); hold on
% errorbar([10:10:40]-1,data_19_Sh5_St63,data_19_Sh5_St63_std,'gs-','linewidth',2); 
% errorbar([10:10:40]+0,data_19_Sh13_St0,data_19_Sh13_St0_std,'ks-','linewidth',2); 
% errorbar([10:10:40]+1,data_19_Sh13_St63,data_19_Sh13_St63_std,'bs-','linewidth',2); 
% 
% 
% set(gca, 'XTick', [10 20 30 40], ...
%          'XTickLabel', [],...
%          'fontSize',16)
% xlim([5 45])
% ylim([0 4])
% % xlabel('bottle size [mL]')
% if iLoc ==1; ylabel('N_2O [\mug N/L]'); end
% title([location{iLoc}])
% 
% 
% if iLoc == 3
%     legend(lgndTxt{1},lgndTxt{2},lgndTxt{3},lgndTxt{4},'interpreter','none','location','north')
%     legend box off
% end
% if iLoc ~= 1
%     set(gca,'YTickLabel',[])
% end
% 
% 
% 
% 
% subplot(6,18,locMat{iLoc,2}(:))
% errorbar([10:10:40]-2,data_441_Sh5_St0,data_441_Sh5_St0_std,'ro--','linewidth',2); hold on
% errorbar([10:10:40]-1,data_441_Sh5_St63,data_441_Sh5_St63_std,'go--','linewidth',2); 
% errorbar([10:10:40]+0,data_441_Sh13_St0,data_441_Sh13_St0_std,'ko--','linewidth',2); 
% errorbar([10:10:40]+1,data_441_Sh13_St63, data_441_Sh13_St63_std,'bo--','linewidth',2); 
% 
% set(gca, 'XTick', [10 20 30 40], ...
%          'XTickLabel', bottleSize,...
%          'fontSize',16)
% xlim([5 45])
% ylim([0 4])
% xlabel('bottle size [mL]')
% if iLoc ==1; ylabel('N_2O [\mug N/L]'); end
% % title([location{iLoc}])
% 
% 
% if iLoc == 3
%     legend(lgndTxt{5},lgndTxt{6},lgndTxt{7},lgndTxt{8},'interpreter','none','location','north')
%     legend box off
% end
% if iLoc ~= 1
%     set(gca,'YTickLabel',[])
% end
% 
% 
% % set(gca,'XTicks', 1:4)
% 
% end
% 
% 
% fig = gcf;
% fig.Units = 'inches';
% fig.PaperPosition = [0 0 24 8];
% 
% flNameFig = sprintf('%s/processed_hsExtract/Apr2026/fig_bottleSize_err.png',pwd);
% 
% 
% print(fig,flNameFig,'-dpng')
% 
% % boxplot([data_19_Sh5_St0(:) data_441_Sh5_St0(:)])% dataW3_19(:) dataW3_441(:) dataW4_441(:) dataW4_19(:)])
% % hold on;
% %%
% % plot([0 7],[mean(n2o_atm_avg_vol_time(:)) mean(n2o_atm_avg_vol_time(:))])
% % 
% % yl = ylim;
% % % text(1,mean(n2o_atm_avg_vol_time(:))+0.02*diff(yl),'fontsize',14)
% % text([1:6]-0.4,mean([dataRR_19(:) dataRR_441(:) dataW3_19(:) dataW3_441(:) dataW4_441(:) dataW4_19(:)]-0.01*diff(yl)),...
% %     {'N_2','He','N_2','He','N_2','He',},'fontsize',14)
% % 
% % 
% % ylabel('N_2O [ug N/L]')
% % ax = gca;
% % 
% % ax.XTick = [1.5 3.5 5.5]
% % ax.XTickLabels = {'Rise', 'Well 3', 'Well 4'}
% % ax.FontSize = 14;
% % 
% % 
% % fig = gcf;
% % fig.Units = 'inches';
% % fig.PaperPosition = [0 0 7 8];
% % flNameFig = sprintf('He_N2_comparision.png');
% % print(fig,flNameFig,'-dpng')
% 
% 
% %%
% RR_441 = n2o_avg_loc_vol_ratio_sht_stt(1,:,2,:,:);
% RR_19 =  n2o_avg_loc_vol_ratio_sht_stt(1,:,1,:,:);
% 
% W3_441 = n2o_avg_loc_vol_ratio_sht_stt(2,:,2,:,:);
% W3_19 =  n2o_avg_loc_vol_ratio_sht_stt(2,:,1,:,:);
% 
% W4_441 = n2o_avg_loc_vol_ratio_sht_stt(3,:,2,:,:);
% W4_19 =  n2o_avg_loc_vol_ratio_sht_stt(3,:,1,:,:);
% 
% 
% 
% 
% [ mean(RR_441(:)) mean(RR_19(:)) mean(RR_441(:))/mean(RR_19(:)) ]
% [mean(W3_441(:)) mean(W3_19(:)) mean(W3_441(:))/mean(W3_19(:))]
% [mean(W4_441(:)) mean(W4_19(:)) mean(W4_441(:))/mean(W4_19(:))]
%%%% plot mean N2O for location * volume * shaking time
% clf; clc
% 
% [xx,yy] = meshgrid(1:3, 1:4);
% az = -80;
% el = 15;
% 
% for iGas = 1:2
%     clf
% 
%     fig = gcf;
%     fig.Units = 'inches';
%     fig.PaperPosition = [0 0 7 8];
% 
%     flNameFig = sprintf('HS_extract_%s.png', inertGas{iGas});
% 
%     %-----------------------------
%     % Data
%     %-----------------------------
%     data = squeeze(n2o_avg_loc_vol_ratio_sht_stt(iGas,1:3,:,:));   % 3 x 4 x 3
%     alphaVals = [1 0.8 0.8];
% 
%     % Center each plane by its own mean for coloring
%     cdata = zeros(size(data));
%     for k = 1:3
%         cdata(k,:,:) = data(k,:,:) - mean(data(k,:,:), 'all');
%     end
% 
%     clim_common = max(abs(cdata(:))) * [-1 1];
% 
%     %-----------------------------
%     % Plot
%     %-----------------------------
%     ax = axes;
%     hold(ax,'on')
% 
%     for k = 1:3
%         surf(ax, xx, yy, squeeze(data(k,:,:)), squeeze(cdata(k,:,:)), ...
%             'FaceColor','interp', ...
%             'EdgeColor','none', ...
%             'FaceAlpha',alphaVals(k));
%     end
% 
%     view(ax, az, el)
%     colormap(ax, jet)
%     clim(ax, clim_common)
% 
%     xlim(ax,[0.7 3.3])
%     ylim(ax,[0.7 4.3])
% 
%     ax.XTick = 1:3;
%     ax.XTickLabel = string(shakeTime);
%     ax.YTick = 1:4;
%     ax.YTickLabel = string(bottleSize);
%     ax.FontSize = 14;
%     ax.Position = [0.12 0.12 0.62 0.78];
% 
%     ylabel(ax,'Volume [ml]')
%     zlabel(ax,'N_2O [\mug N/L]')
%     title(ax, sprintf('Inert Gas: %s', inertGas{iGas}), 'FontSize', 16)
% 
%     box(ax,'on')
%     grid(ax,'on')
% 
%     %-----------------------------
%     % Colorbar
%     %-----------------------------
%     cbh = colorbar(ax, 'eastoutside');
%     cbh.Position = [0.85 0.20 0.025 0.70];
%     cbh.Ticks = linspace(clim_common(1), clim_common(2), 7);
%     cbh.TickLabels = compose('%.3f', cbh.Ticks);
% 
%     %-----------------------------
%     % Custom x-label as rotated text
%     %-----------------------------
%     drawnow
%     xl = xlim(ax);
%     yl = ylim(ax);
%     zl = zlim(ax);
% 
%     text(ax, mean(xl), yl(1) - 0.10*range(yl), zl(1), 'Shake time [min]', ...
%         'HorizontalAlignment','center', ...
%         'VerticalAlignment','top', ...
%         'FontSize',14, ...
%         'Rotation', -(az + 13));
% 
%     print(fig, flNameFig, '-dpng')
% end
