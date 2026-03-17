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
flNameSample = sprintf('%s/processed_hsExtract/sample.mat',pwd);
load(flNameSample,'T_sample')

inertGas = {'N2','HE'};

location = {'RR','W3','W4'};
bottleSize = [10 25 72 158]; %[mL]
shakeTime = [3 5 8]; %[min]

%% ======================================================================== 
tempK = T_sample.tempK;
peakH = T_sample.peakH;

bottleVol = T_sample.bottleSize/1000; %[L]
Vaq = bottleVol*9/10; %[L] 
Vg =  bottleVol*1/10; %[L] 
Vsyringe = 1/1000; %[L]

valveVol = 3/1000; %[L]
n2o_ppm_valve = 0/1000; %[ppm]

logK0 = A + (B./tempK) + C*(log(tempK));
K0 = exp(logK0); % Henry's law constant [mol/L/atm]
    
    % Not valve corrected -----------------------------------------
    n2o_ppm_uncorrected = (peakH*slope + offset);
    
    % Valve corrected -----------------------------------------
    n2o_ppm = (n2o_ppm_uncorrected .* (Vsyringe+valveVol) )./Vsyringe;
    
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
    logK0_atm = A + (B./mean(tempK)) + C*(log(mean(tempK)));
    K0_atm = exp(logK0_atm); % Henry's law constant [mol/L/atm]
    
    n2o_atm_atmUnit = n2o_ppm_valve * 1e-6;
    ng_atm = (n2o_atm_atmUnit .* Vg)./(R*mean(tempK)); % N2O moles in gas phase (from ideal gas law)
    naq_atm = K0_atm.*n2o_atm_atmUnit.*Vaq; % dissolved N2O moles in aqueous phase (from Henry's law)
    
    n2o_atm_mol_L = (ng_atm+naq_atm)./Vaq; % Total N2O concentration in mol/L
    n2o_atm_molN_L = n2o_atm_mol_L * 2;
    n2o_atm_gN_L = n2o_atm_molN_L * amuN;
    n2o_atm_mgN_L = n2o_atm_gN_L * 1e+3;
    n2o_atm_ugN_L = n2o_atm_mgN_L * 1e+3;

% into data table
% T_sample = removevars(T_sample,'n2o_ppm_measured');
T_sample.n2o_ppm_uncorrected = n2o_ppm_uncorrected;
T_sample.n2o_ppm = n2o_ppm;
T_sample.n2o_ugN_L = n2o_ugN_L;
T_sample.n2o_atm_ugN_L = n2o_atm_ugN_L;



%% data arrangements

n2o_avg_vol_time_uncorrected = size(numel(inertGas),numel(location),numel(bottleSize),numel(shakeTime));
n2o_std_vol_time_uncorrected = size(numel(inertGas),numel(location),numel(bottleSize),numel(shakeTime));

n2o_avg_vol_time = size(numel(inertGas),numel(location),numel(bottleSize),numel(shakeTime));
n2o_std_vol_time = size(numel(inertGas),numel(location),numel(bottleSize),numel(shakeTime));

n2o_atm_avg_vol_time = size(numel(inertGas),numel(location),numel(bottleSize),numel(shakeTime));

for iGas = 1:numel(inertGas)
for iLoc = 1:numel(location)
for iVol = 1:numel(bottleSize)
for iTime = 1:numel(shakeTime)
    tmp = T_sample(strcmp(T_sample.inertGas,inertGas{iGas})... 
           & strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ...
           ,:).n2o_ugN_L;

   

    n2o_avg_vol_time(iGas,iLoc,iVol,iTime) = mean(tmp);
    n2o_std_vol_time(iGas,iLoc,iVol,iTime) = std(tmp);

    n2o_atm_avg_vol_time(iGas,iLoc,iVol,iTime) = mean(n2o_atm_ugN_L( strcmp(T_sample.inertGas,inertGas{iGas})... 
           & strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ));
    
    
end
end
end
end

writetable(T_sample,'HS_extract.xlsx')

%% plot mean N2O for location * volume * shaking time
clf; clc

[xx,yy] = meshgrid(1:3, 1:4);
az = -80;
el = 15;

for iGas = 1:2
    clf

    fig = gcf;
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 7 8];

    flNameFig = sprintf('HS_extract_%s.png', inertGas{iGas});

    %-----------------------------
    % Data
    %-----------------------------
    data = squeeze(n2o_avg_vol_time(iGas,1:3,:,:));   % 3 x 4 x 3
    alphaVals = [1 0.8 0.8];

    % Center each plane by its own mean for coloring
    cdata = zeros(size(data));
    for k = 1:3
        cdata(k,:,:) = data(k,:,:) - mean(data(k,:,:), 'all');
    end

    clim_common = max(abs(cdata(:))) * [-1 1];

    %-----------------------------
    % Plot
    %-----------------------------
    ax = axes;
    hold(ax,'on')

    for k = 1:3
        surf(ax, xx, yy, squeeze(data(k,:,:)), squeeze(cdata(k,:,:)), ...
            'FaceColor','interp', ...
            'EdgeColor','none', ...
            'FaceAlpha',alphaVals(k));
    end

    view(ax, az, el)
    colormap(ax, jet)
    clim(ax, clim_common)

    xlim(ax,[0.7 3.3])
    ylim(ax,[0.7 4.3])

    ax.XTick = 1:3;
    ax.XTickLabel = string(shakeTime);
    ax.YTick = 1:4;
    ax.YTickLabel = string(bottleSize);
    ax.FontSize = 14;
    ax.Position = [0.12 0.12 0.62 0.78];

    ylabel(ax,'Volume [ml]')
    zlabel(ax,'N_2O [\mug N/L]')
    title(ax, sprintf('Inert Gas: %s', inertGas{iGas}), 'FontSize', 16)

    box(ax,'on')
    grid(ax,'on')

    %-----------------------------
    % Colorbar
    %-----------------------------
    cbh = colorbar(ax, 'eastoutside');
    cbh.Position = [0.85 0.20 0.025 0.70];
    cbh.Ticks = linspace(clim_common(1), clim_common(2), 7);
    cbh.TickLabels = compose('%.3f', cbh.Ticks);

    %-----------------------------
    % Custom x-label as rotated text
    %-----------------------------
    drawnow
    xl = xlim(ax);
    yl = ylim(ax);
    zl = zlim(ax);

    text(ax, mean(xl), yl(1) - 0.10*range(yl), zl(1), 'Shake time [min]', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',14, ...
        'Rotation', -(az + 13));

    print(fig, flNameFig, '-dpng')
end

%%
% [inertGas location bottleSize shakeTime]
indBottleSize = 1:4; %[10 25 72 158]
iShakeTime = 1:3; %[3 5 8]
clf

dataRR_N2 = squeeze(n2o_avg_vol_time(1,1,indBottleSize,iShakeTime));
dataRR_HE = squeeze(n2o_avg_vol_time(2,1,indBottleSize,iShakeTime));
dataW3_N2 = squeeze(n2o_avg_vol_time(1,2,indBottleSize,iShakeTime));
dataW3_HE = squeeze(n2o_avg_vol_time(2,2,indBottleSize,iShakeTime));
dataW4_N2 = squeeze(n2o_avg_vol_time(1,3,indBottleSize,iShakeTime));
dataW4_HE = squeeze(n2o_avg_vol_time(2,3,indBottleSize,iShakeTime));

boxplot([dataRR_N2(:) dataRR_HE(:) dataW3_N2(:) dataW3_HE(:) dataW4_HE(:) dataW4_N2(:)])
hold on;

plot([0 7],[mean(n2o_atm_avg_vol_time(:)) mean(n2o_atm_avg_vol_time(:))])

yl = ylim;
% text(1,mean(n2o_atm_avg_vol_time(:))+0.02*diff(yl),'fontsize',14)
text([1:6]-0.4,mean([dataRR_N2(:) dataRR_HE(:) dataW3_N2(:) dataW3_HE(:) dataW4_HE(:) dataW4_N2(:)]-0.01*diff(yl)),...
    {'N_2','He','N_2','He','N_2','He',},'fontsize',14)

                        
ylabel('N_2O [ug N/L]')
ax = gca;

ax.XTick = [1.5 3.5 5.5]
ax.XTickLabels = {'Rise', 'Well 3', 'Well 4'}
ax.FontSize = 14;


fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 8];
flNameFig = sprintf('He_N2_comparision.png');
print(fig,flNameFig,'-dpng')
