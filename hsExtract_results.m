clear all; clc; clf

%% He & N2
% flNameSample = sprintf('%s/processed_HE_N2/HE_N2.mat',pwd);
% load(flNameSample,'T_HE_N2')
% 
% T_HE = T_HE_N2(strcmp(T_HE_N2.gasType, 'HE'),:);
% T_N2 = T_HE_N2(strcmp(T_HE_N2.gasType, 'N2'),:);
% 
% HE_N2O = T_HE.n2o_ugN_L(2:end);
% N2_N2O = T_N2.n2o_ugN_L(2:end);
% 
% M = [mean(HE_N2O) mean(N2_N2O)];
% E = [std(HE_N2O) std(N2_N2O)];
% 
% clf
% errorbar(M,E,'x')
% 
% xlim([0.5 2.5])
% 
% ax = gca;
% ax.XTick = [1 2];
% ax.XTickLabel = {"He","N2"}
% ax.FontSize = 20;
% ylabel('N2O [ug/ml]')
% xlabel('Gas type')


%% load sample data
flNameSample = sprintf('%s/processed_hsExtract/sample.mat',pwd);
load(flNameSample,'T_sample')

T_rise = T_sample(strcmp(T_sample.location, 'RISE'),:);

location = {'RISE','WELL3','WELL4'};
bottleSize = [7.84 24 70 157];
shakeTime = [3 5 8];

n2o_avg_vol_time = size(numel(location),numel(bottleSize),numel(shakeTime));
n2o_std_vol_time = size(numel(location),numel(bottleSize),numel(shakeTime));

for iLoc = 1:numel(location)
for iVol = 1:numel(bottleSize)
for iTime = 1:numel(shakeTime)
    tmp = T_sample(strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ...
           & ~T_sample.suspicious...
           ,:).n2o_ugN_L;

    % tmp = T_sample(strcmp(T_sample.location,location{iLoc}) ...
    %        & T_sample.bottleSize==bottleSize(iVol) ...
    %        & T_sample.shakeTime==shakeTime(iTime) ...
    %        ,:).n2o_ugN_L;

    n2o_avg_vol_time(iLoc,iVol,iTime) = mean(tmp);
    n2o_std_vol_time(iLoc,iVol,iTime) = std(tmp);
    
end
end
end

%%
flNameSTD = sprintf('%s/processed/STD.mat',pwd);
load(flNameSTD,'T_STD','slope','offset')
%% ========================================================================
clc
A = -148.1; B = 8610; C = 20.266; 
R = 0.0821; %[L*atm/K/mol]
amuN = 14.0067; 
amuO = 16;
amuN2O = 2*amuN + amuO;
% 
tempK = T_sample.tempK;
peakH = T_sample.peakH;
dilutionFactor = T_sample.dilutionFactor;
bottleSize = T_sample.bottleSize;

bottleVol = bottleSize/1000; %[L]
Vaq = bottleVol*9/10; %[L] 
Vg =  bottleVol*1/10; %[L] 

valveVol = 3/1000; %[L]
n2o_ppm_atm = 340/1000; %[ppm]

logK0 = A + (B./tempK) + C*(log(tempK));
K0 = exp(logK0); % Henry's law constant [mol/L/atm]
    
    % Not valve corrected -----------------------------------------
    n2o_ppm_uncorrected = (peakH*slope + offset).*dilutionFactor;
    n2o_atmUnit_uncorrected = n2o_ppm_uncorrected * 1e-6;

    ng_uncorrected = (n2o_atmUnit_uncorrected .* Vg)./(R*tempK); % N2O moles in gas phase (from ideal gas law)
    naq_uncorrected = K0.*n2o_atmUnit_uncorrected.*Vaq; % dissolved N2O moles in aqueous phase (from Henry's law)
    
    n2o_mol_L_uncorrected = (ng_uncorrected+naq_uncorrected)./Vaq; % Total N2O concentration in mol/L
    n2o_molN_L_uncorrected = n2o_mol_L_uncorrected * 2;
    n2o_gN_L_uncorrected = n2o_molN_L_uncorrected * amuN;
    n2o_mgN_L_uncorrected = n2o_gN_L_uncorrected * 1e+3;
    n2o_ugN_L_uncorrected = n2o_mgN_L_uncorrected * 1e+3;

    % Valve corrected -----------------------------------------
    n2o_ppm = (n2o_ppm_uncorrected .* (bottleSize/1000+valveVol) - (valveVol*n2o_ppm_atm))./bottleVol;
    
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
    
    n2o_atm_atmUnit = n2o_ppm_atm * 1e-6;
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
T_sample.n2o_molN_L_uncorrected = n2o_molN_L_uncorrected;
T_sample.n2o_ugN_L_uncorrected = n2o_ugN_L_uncorrected;
T_sample.n2o_ugN_L = n2o_ugN_L;
T_sample.n2o_atm_ugN_L = n2o_atm_ugN_L;


bottleSize = [7.84 24 70 157];
shakeTime = [3 5 8];

n2o_avg_vol_time_uncorrected = size(numel(location),numel(bottleSize),numel(shakeTime));
n2o_std_vol_time_uncorrected = size(numel(location),numel(bottleSize),numel(shakeTime));

n2o_avg_vol_time = size(numel(location),numel(bottleSize),numel(shakeTime));
n2o_std_vol_time = size(numel(location),numel(bottleSize),numel(shakeTime));

n2o_atm_avg_vol_time = size(numel(location),numel(bottleSize),numel(shakeTime));


for iLoc = 1:numel(location)
for iVol = 1:numel(bottleSize)
for iTime = 1:numel(shakeTime)
    tmp = T_sample(strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ...
           & ~T_sample.suspicious...
           ,:).n2o_ugN_L;

    tmp_uncorrected = T_sample(strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ...
           & ~T_sample.suspicious...
           ,:).n2o_ugN_L_uncorrected;

    % tmp = T_sample(strcmp(T_sample.location,location{iLoc}) ...
    %        & T_sample.bottleSize==bottleSize(iVol) ...
    %        & T_sample.shakeTime==shakeTime(iTime) ...
    %        ,:).n2o_ugN_L;

    n2o_avg_vol_time(iLoc,iVol,iTime) = mean(tmp);
    n2o_std_vol_time(iLoc,iVol,iTime) = std(tmp);

    n2o_avg_vol_time_uncorrected(iLoc,iVol,iTime) = mean(tmp_uncorrected);
    n2o_std_vol_time_uncorrected(iLoc,iVol,iTime) = std(tmp_uncorrected);

    n2o_atm_avg_vol_time(iLoc,iVol,iTime) = mean(n2o_atm_ugN_L(strcmp(T_sample.location,location{iLoc}) ...
           & T_sample.bottleSize==bottleSize(iVol) ...
           & T_sample.shakeTime==shakeTime(iTime) ...
           & ~T_sample.suspicious));
    
    
end
end
end


writetable(T_sample,'HS_extract.xlsx')

%% valve corrected N2O
clf; clc

flNameFig = sprintf('HS_extract_corrected.png');

% figure
clr = 'rkc';
% for iLoc = 1:numel(location)

[xx,yy] = meshgrid(1:3, 1:4);

fx1= -50; 
fx2= 13.58;

% subplot(121)
data1 = squeeze(n2o_avg_vol_time(1,:,:));
ax1 = axes;
surf(ax1,xx,yy,squeeze(data1), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 1);
cl1 = clim(ax1);
mean1 = mean(cl1);

view(fx1,fx2)
cbh = colorbar
cbh.Ticks = linspace(cl1(1),cl1(2),7);
cbh.TickLabels = sprintf('%.3f\n',cbh.Ticks' - mean1);


data2 = squeeze(n2o_avg_vol_time(2,:,:));
ax2 = axes;
surf(ax2,xx,yy,squeeze(data2), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
cl2 = clim(ax2);
mean2 = mean(cl2);
view(fx1,fx2)
clim(ax2,cl1-mean1+mean2)
axis off


data3 = squeeze(n2o_avg_vol_time(3,:,:));
ax3 = axes;
surf(ax3,xx,yy,data3, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
cl3 = clim(ax3);
mean3 = mean(cl3);
view(fx1,fx2)
clim(ax3,cl1-mean1+mean3)
axis off


data4 = squeeze(n2o_atm_avg_vol_time(3,:,:));
ax4 = axes;
surf(ax4,xx,yy,data4, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
cl4 = clim(ax4);
mean4 = mean(cl4);
view(fx1,fx2)
clim(ax4,cl1-mean1+mean4)
axis off


ax2.Position = ax1.Position;
ax3.Position = ax1.Position;
ax4.Position = ax1.Position;

ax2.Color = 'none'; % set the background color for ax2 to transparent
ax3.Color = 'none'; % set the background color for ax2 to transparent
ax4.Color = 'none'; % set the background color for ax2 to transparent

colmap = 'jet';
colormap(ax1,colmap); %colorbar; % change the color of s1
colormap(ax2,colmap); % change the color of s2
colormap(ax3,colmap); % change the color of s2
colormap(ax4,colmap); % change the color of s2

linkaxes([ax1 ax2 ax3 ax4])

xlabel('shake time [min]'); xlim([0.7 3.3])
ylabel('volume [ml]'); ylim([0.7 4.3])
zlabel('N_2O [ug N/L]'); zlim([0.2 2])

ax1.XTick = 1:3; ax1.XTickLabel = num2cell(shakeTime);
ax1.YTick = 1:4; ax1.YTickLabel = num2cell(bottleSize);

xlabel(ax1,'Shake time [min]')
ylabel(ax1,'Volume [ml]')
zlabel(ax1,'N_2O [μg N/l]')
set(ax1,'fontSize',14)

fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 8];

print(fig,flNameFig,'-dpng')

% shading interp
% % shading flat

'done'

%% both
% clf; clc
% 
% flNameFig = sprintf('HS_extract_both.png');
% 
% % figure
% clr = 'rkc';
% % for iLoc = 1:numel(location)
% 
% [xx,yy] = meshgrid(1:3, 1:4);
% 
% fx1= -50; 
% fx2= 13.58;
% 
% % subplot(121)
% data11 = squeeze(n2o_avg_vol_time_uncorrected(1,:,:));
% ax11 = axes;
% surf(ax11,xx,yy,squeeze(data11), 'FaceColor', 'interp', 'EdgeColor','none',...
%                                  'FaceAlpha', 0.7);
% cl11 = clim(ax11);
% mean11 = mean(cl11);
% 
% view(fx1,fx2)
% cbh = colorbar;
% cbh.Ticks = linspace(cl11(1),cl11(2),7);
% cbh.TickLabels = sprintf('%.3f\n',cbh.Ticks' - mean11);
% 
% 
% data1 = squeeze(n2o_avg_vol_time(1,:,:));
% ax1 = axes;
% surf(ax1,xx,yy,squeeze(data1), 'FaceColor', 'interp', 'EdgeColor','none',...
%                                'FaceAlpha', 0.7);
% cl1 = clim(ax1);
% mean1 = mean(cl1);
% view(fx1,fx2)
% clim(ax1,cl11-mean11+mean1)
% axis off
% 
% 
% data22 = squeeze(n2o_avg_vol_time_uncorrected(2,:,:));
% ax22 = axes;
% surf(ax22,xx,yy,squeeze(data22), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% cl22 = clim(ax22);
% mean22 = mean(cl22);
% view(fx1,fx2)
% clim(ax22,cl11-mean11+mean22)
% axis off
% 
% 
% data2 = squeeze(n2o_avg_vol_time(2,:,:));
% ax2 = axes;
% surf(ax2,xx,yy,squeeze(data2), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% cl2 = clim(ax2);
% mean2 = mean(cl2);
% view(fx1,fx2)
% clim(ax2,cl11-mean11+mean2)
% axis off
% 
% 
% data33 = squeeze(n2o_avg_vol_time_uncorrected(3,:,:));
% ax33 = axes;
% surf(ax33,xx,yy,squeeze(data33), 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% cl33 = clim(ax33);
% mean33 = mean(cl33);
% view(fx1,fx2)
% clim(ax33,cl11-mean11+mean33)
% axis off
% 
% 
% 
% data3 = squeeze(n2o_avg_vol_time(3,:,:));
% ax3 = axes;
% surf(ax3,xx,yy,data3, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
% cl3 = clim(ax3);
% mean3 = mean(cl3);
% view(fx1,fx2)
% clim(ax3,cl11-mean11+mean3)
% axis off
% 
% 
% data4 = squeeze(n2o_atm_avg_vol_time(3,:,:));
% ax4 = axes;
% surf(ax4,xx,yy,data4, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
% cl4 = clim(ax4);
% mean4 = mean(cl4);
% view(fx1,fx2)
% clim(ax4,cl11-mean11+mean4)
% axis off
% 
% ax1.Position = ax11.Position;
% ax22.Position = ax11.Position;
% ax2.Position = ax11.Position;
% ax33.Position = ax11.Position;
% ax3.Position = ax11.Position;
% ax4.Position = ax11.Position;
% 
% ax1.Color = 'none'; % set the background color for ax2 to transparent
% ax2.Color = 'none'; % set the background color for ax2 to transparent
% ax3.Color = 'none'; % set the background color for ax2 to transparent
% ax4.Color = 'none'; % set the background color for ax2 to transparent
% 
% colmap = 'jet';
% colormap(ax11,colmap); %colorbar; % change the color of s1
% colormap(ax1,colmap); %colorbar; % change the color of s1
% colormap(ax2,colmap); % change the color of s2
% colormap(ax22,colmap); % change the color of s2colormap(ax2,colmap); % change the color of s2
% colormap(ax3,colmap); % change the color of s2
% colormap(ax33,colmap); % change the color of s2
% colormap(ax4,colmap); % change the color of s2
% 
% linkaxes([ax11 ax22 ax33 ax1 ax2 ax3 ax4])
% 
% xlabel('shake time [min]'); xlim([0.7 3.3])
% ylabel('volume [ml]'); ylim([0.7 4.3])
% zlabel('N_2O [ug N/L]'); zlim([0.2 2])
% 
% ax11.XTick = 1:3; ax11.XTickLabel = num2cell(shakeTime);
% ax11.YTick = 1:4; ax11.YTickLabel = num2cell(bottleSize);
% 
% xlabel(ax11,'Shake time [min]')
% ylabel(ax11,'Volume [ml]')
% zlabel(ax11,'N_2O [μg N/l]')
% set(ax11,'fontSize',14)
% 
% fig = gcf;
% fig.Units = 'inches';
% fig.PaperPosition = [0 0 7 8];
% 
% print(fig,flNameFig,'-dpng')
% 
% % shading interp
% % % shading flat
% 
% 'done'

%%
clf

data11 = squeeze(n2o_avg_vol_time_uncorrected(1,:,:));
data1 = squeeze(n2o_avg_vol_time(1,:,:));
data22 = squeeze(n2o_avg_vol_time_uncorrected(2,:,:));
data2 = squeeze(n2o_avg_vol_time(2,:,:));
data33 = squeeze(n2o_avg_vol_time_uncorrected(3,:,:));
data3 = squeeze(n2o_avg_vol_time(3,:,:));
data4 = squeeze(n2o_atm_avg_vol_time(3,:,:));

boxplot([data11(:) data1(:) data22(:) data2(:) data3(:) data33(:)])
hold on;

plot([0 7],[mean(n2o_atm_avg_vol_time(:)) mean(n2o_atm_avg_vol_time(:))])

yl = ylim;
text(1,mean(n2o_atm_avg_vol_time(:))+0.02*diff(yl),'ATM','fontsize',14)
text([1:6]-0.5,min([data11(:) data1(:) data22(:) data2(:) data3(:) data33(:)]-0.01*diff(yl)),...
    {'uncorrected','valve-corrected','uncorrected','valve-corrected','uncorrected','valve-corrected'},...
    'fontsize',10)

                        
ylabel('N_2O [ug N/L]')
ax = gca;

ax.XTick = [1.5 3.5 5.5]
ax.XTickLabels = {'Rise', 'Well 3', 'Well 4'}
ax.FontSize = 14;


fig = gcf;
fig.Units = 'inches';
    fig.PaperPosition = [0 0 7 8];
    flNameFig = sprintf('HS_extract_box.png');
    print(fig,flNameFig,'-dpng')
