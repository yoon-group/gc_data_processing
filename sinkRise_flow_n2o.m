clf; clear all; clc
%% ========================================================================
% 
%                   2024 O'leno park fieldwork 
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% This code is to visualize the flow conditions, SW-GW interaction, and N2O 
% production at O'leno state park in 2024 Fall. This code download the flow 
% rates at the Sink and Rise USGS gauges from the USGS data server.
% =========================================================================

%% download USGS data 
startDT = "2025-1-1 00:00"; endDT = "2025-12-31 00:00";
% startDT = "2025-1-15 00:00"; endDT = "2025-3-2 00:00";
% 
parameterCd = "00060"; %streamflow cubic feet per second
createFig = false;

a.sites = "02321898"; a.startDT = startDT; a.endDT = endDT; 
a.parameterCd = parameterCd; a.createFig = createFig;
[sink] = scrape(a);


a.sites = "02321958"; a.startDT = startDT; a.endDT = endDT; 
a.parameterCd = parameterCd; a.createFig = createFig ;
[rise] = scrape(a);


%% dye tracing experiment time
t_tracingStart = datetime('2024-11-18T14:10','InputFormat','uuuu-MM-dd''T''HH:mm','TimeZone','local');
t_tracingEnd = datetime('2024-12-02T12:10','InputFormat','uuuu-MM-dd''T''HH:mm','TimeZone','local');


%% sampling dates & temperature data
% ================= sampled from the shore near the steps =================
% dateFieldWork(1,1) = datetime('10-04-2024 10:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
% dateFieldWork(end+1,1) = datetime('10-14-2024 13:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
% =========================================================================
dateFieldWork(1,1) = datetime('10-16-2024 13:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('10-17-2024 10:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local');
dateFieldWork(end+1,1) = datetime('10-21-2024 10:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local');
dateFieldWork(end+1,1) = datetime('10-24-2024 09:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local');
dateFieldWork(end+1,1) = datetime('11-14-2024 09:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local');
% dateFieldWork(end+1,1) = datetime('12-2-2024 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); % inconsistent
dateFieldWork(end+1,1) = datetime('12-17-2024 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
% dateFieldWork(end+1,1) = datetime('01-02-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('01-24-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('01-27-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('01-30-2025 14:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('02-03-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('02-10-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 
dateFieldWork(end+1,1) = datetime('02-21-2025 11:30:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local'); 

dateFieldWork(find(dateFieldWork<rise.dte(1))) = [];
%% load sample data
clc
flNameSample = sprintf('%s/processed/sample.mat',pwd);
load(flNameSample,'T_sample')
T_n2o = T_sample; clear T_sample

flNameYSI_NO3 = sprintf('%s/N2O_master.xlsx',pwd);
T_GW = readtable(flNameYSI_NO3,'Sheet','OlenoWells');

T_SW = readtable(flNameYSI_NO3,'Sheet','Sweet-Rise trends');


flName_dPCR = sprintf('/home/public/microbiomeData/N2O_field_sampling_Feb_2025.csv');
T_dPCR = readtable(flName_dPCR);

flName16S_rDNA = sprintf('/home/public/microbiomeData/wf_16s_rDNA_counts_genus_02212025.xlsx');
T_16S_rDNA_well3 = readtable(flName16S_rDNA,'Sheet','WELL3',ReadVariableNames=true);
T_16S_rDNA_well4 = readtable(flName16S_rDNA,'Sheet','WELL4',ReadVariableNames=true);

%% exact time matching in T_16S_rDNA
varNames = T_16S_rDNA_well3.Properties.VariableNames(2:end);
for iDate = 1:length(varNames)
    txt = datetime(replace(varNames{iDate}(2:end),'_','-'));
    tmp = datetime(txt);
    tmp.TimeZone = 'local';
    date_16s_rDNA_well3(iDate,1) = tmp;
end
for iRow = 1:size(date_16s_rDNA_well3,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           date_16s_rDNA_well3(iRow) == 0);
    if indField
        date_16s_rDNA_well3(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end


varNames = T_16S_rDNA_well4.Properties.VariableNames(2:end);
for iDate = 1:length(varNames)
    txt = datetime(replace(varNames{iDate}(2:end),'_','-'));
    tmp = datetime(txt);
    tmp.TimeZone = 'local';
    date_16s_rDNA_well4(iDate,1) = tmp;
end
for iRow = 1:size(date_16s_rDNA_well4,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           date_16s_rDNA_well4(iRow) == 0);
    if indField
        date_16s_rDNA_well4(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end

%%
clc
for location = {'well3','well4'}
    eval(sprintf('microbeName  = T_16S_rDNA_%s{:,1};',location{:}));
    eval(sprintf('T_16S_rDNA = T_16S_rDNA_%s{:,2:end};',location{:}));
        
    totalCount = sum(T_16S_rDNA);
    portionRate = T_16S_rDNA./totalCount;
    portionRate_avg = mean(portionRate,2);
    
    [~,ind] = sort(portionRate_avg,'descend');
    
    eval(sprintf('microbeName_%s = microbeName(ind,:);',location{:}));
    eval(sprintf('portionRate_%s = portionRate(ind,:);',location{:}));
    
end

%% exact time matching in T_dPCR
clear dPcrDate
for iRow = 1:size(T_dPCR,1)
    tmp = T_dPCR.Date(iRow);
    tmp.TimeZone = 'local';
    dPcrDate(iRow,1) = tmp;
end
T_dPCR.Date = dPcrDate;

for iRow = 1:size(T_dPCR,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           T_dPCR.Date(iRow) == 0);
    if indField
        T_dPCR.sampleDate(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end

%% exact time matching in T_n2o
for iRow = 1:size(T_n2o,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           T_n2o.sampleDate(iRow) == 0);
    if indField
        T_n2o.sampleDate(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end

%% exact time matching in T_GW
clear ysiDate
for iRow = 1:size(T_GW,1)
    tmp = T_GW.Date(iRow);
    tmp.TimeZone = 'local';
    ysiDate(iRow,1) = tmp;
end
T_GW.Date = ysiDate;

for iRow = 1:size(T_GW,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           T_GW.Date(iRow) == 0);
    if indField
        T_GW.sampleDate(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end

%% exact time matching in T_SW
clear ysiDate
for iRow = 1:size(T_SW,1)
    tmp = T_SW.SampleDate(iRow);
    tmp.TimeZone = 'local';
    ysiDate(iRow,1) = tmp;
end
T_SW.Date = ysiDate;

for iRow = 1:size(T_SW,1)
    indField = find(dateFieldWork-timeofday(dateFieldWork)- ...
                           T_SW.Date(iRow) == 0);
    if indField
        T_SW.sampleDate(iRow,1) = dateFieldWork(indField);
    end
    clear indField
end


%% YSI / IC / FDOM / CDOM
clc
T_no3 = table(dateFieldWork);
T_DO = table(dateFieldWork);
T_ORP = table(dateFieldWork);
T_SpC = table(dateFieldWork);
T_DOC = table(dateFieldWork);
T_d18O = table(dateFieldWork);
T_dD = table(dateFieldWork);
T_nh4 = table(dateFieldWork);
T_fdom = table(dateFieldWork);
T_fdom_b = table(dateFieldWork);
T_fdom_t = table(dateFieldWork);
T_fdom_a = table(dateFieldWork);
T_fdom_m = table(dateFieldWork);
T_fdom_c = table(dateFieldWork);
T_hix = table(dateFieldWork);

location = {'Well 3','Well 4','Well 4A','Rise','Sweet'};

for iLoc = 1:length(location)

    clear no3_ugN_L  DO_mg_L ORP_mV SpC_uS_cm_ DOC_mgC_L_ ...
            d18O_permil_VsVSMOW_   dD_permil_VsVSMOW_ ...
            NH4_ugN_L_   TotalFDOM_R_U__ ...
            b_R_U_  t_R_U__ a_R_U__  m_R_U__  c_R_U__ ...
            hix

    
    T_data = T_GW;
    if sum(strcmpi({'Rise','Sweet'},location{iLoc}))
        T_data = T_SW;
        T_data = renamevars(T_data,'Site','Location');
        T_data = renamevars(T_data,'NO3__gN_L_','NO3__ugN_L_'); % why not NO3__ugN_L_? -> The leter conversion error, but okay.
        T_data = renamevars(T_data,'SpC__S_cm__YSI_','SpC_uS_cm_'); % why not SpC_uS_cm_? -> The leter conversion error, but okay.
        T_data = renamevars(T_data,'NH4__gN_L_','NH4_ugN_L_');% why not NH4_ugN_L_? -> The leter conversion error, but okay.
    end
    
    for iFieldDate = 1:numel(dateFieldWork)

        ind = find(T_data.sampleDate == dateFieldWork(iFieldDate)...
                   & strcmpi(T_data.Location,location{iLoc}));

        if ~isempty(ind)
            no3_ugN_L(iFieldDate,1) = T_data.NO3__ugN_L_(ind);
            DO_mg_L(iFieldDate,1) = T_data.DO_mg_L_(ind);
            ORP_mV(iFieldDate,1) = T_data.ORP_mV_(ind);
            SpC_uS_cm_(iFieldDate,1) = T_data.SpC_uS_cm_(ind);
            DOC_mgC_L_(iFieldDate,1) = T_data.DOC_mgC_L_(ind);
            d18O_permil_VsVSMOW_(iFieldDate,1) = T_data.d18O_permil_VsVSMOW_(ind);
            dD_permil_VsVSMOW_(iFieldDate,1) = T_data.dD_permil_VsVSMOW_(ind);
            NH4_ugN_L_(iFieldDate,1) = T_data.NH4_ugN_L_(ind);
            b_R_U_(iFieldDate,1) = T_data.b_R_U_(ind);
            t_R_U__(iFieldDate,1) = T_data.t_R_U__(ind);
            a_R_U__(iFieldDate,1) = T_data.a_R_U__(ind);
            m_R_U__(iFieldDate,1) = T_data.m_R_U__(ind);
            c_R_U__(iFieldDate,1) = T_data.c_R_U__(ind);
            TotalFDOM_R_U__(iFieldDate,1) = T_data.TotalFDOM_R_U__(ind);
            HIX(iFieldDate,1) = T_data.HIX(ind);
        else
            no3_ugN_L(iFieldDate,1) = nan;
            DO_mg_L(iFieldDate,1) = nan;
            ORP_mV(iFieldDate,1) = nan;
            SpC_uS_cm_(iFieldDate,1) = nan;
            DOC_mgC_L_(iFieldDate,1) = nan;
            d18O_permil_VsVSMOW_(iFieldDate,1) = nan;
            dD_permil_VsVSMOW_(iFieldDate,1) = nan;
            NH4_ugN_L_(iFieldDate,1) = nan;
            b_R_U_(iFieldDate,1) = nan;
            t_R_U__(iFieldDate,1) = nan;
            a_R_U__(iFieldDate,1) = nan;
            m_R_U__(iFieldDate,1) = nan;
            c_R_U__(iFieldDate,1) = nan;
            TotalFDOM_R_U__(iFieldDate,1) = nan;
            HIX(iFieldDate,1) = nan;
        end

    end
    eval(sprintf('T_no3.%s = no3_ugN_L;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_DO.%s = DO_mg_L;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_ORP.%s = ORP_mV;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_SpC.%s = SpC_uS_cm_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_DOC.%s = DOC_mgC_L_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_d18O.%s = d18O_permil_VsVSMOW_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_dD.%s = dD_permil_VsVSMOW_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_nh4.%s = NH4_ugN_L_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom.%s = TotalFDOM_R_U__;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom_b.%s = b_R_U_;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom_t.%s = t_R_U__;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom_a.%s = a_R_U__;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom_m.%s = m_R_U__;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_fdom_c.%s = c_R_U__;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_hix.%s = HIX;',upper(erase(location{iLoc},' '))))

end



%% dPCR readings
clc
T_16s = table(dateFieldWork);
T_Nir = table(dateFieldWork);
T_Nos = table(dateFieldWork);
T_Nir_Nos = table(dateFieldWork);

location = {'Well 3','Well 4'};

for iLoc = 1:length(location)

    clear x16s  Nir Nos Nir_Nos

    for iFieldDate = 1:numel(dateFieldWork)

        ind = find(T_dPCR.sampleDate == dateFieldWork(iFieldDate)...
                   & strcmpi(T_dPCR.Location,location{iLoc}));

        if ~isempty(ind)
            x16s(iFieldDate,1) = T_dPCR.x16s(ind);
            Nir(iFieldDate,1) = T_dPCR.Nir(ind);
            Nos(iFieldDate,1) = T_dPCR.Nos(ind);
            Nir_Nos(iFieldDate,1) = T_dPCR.Nir(ind)/T_dPCR.Nos(ind);
        else
            x16s(iFieldDate,1) = nan;
            Nir(iFieldDate,1) = nan;
            Nos(iFieldDate,1) = nan;
            Nir_Nos(iFieldDate,1) = nan;            
        end

    end
    eval(sprintf('T_16s.%s = x16s;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_Nir.%s = Nir;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_Nos.%s = Nos;',upper(erase(location{iLoc},' '))))
    eval(sprintf('T_Nir_Nos.%s = Nir_Nos;',upper(erase(location{iLoc},' '))))
end
%% mean n2o concentration - ATM [ug N/L]
T_n2o_ugN_L_dATM = table(dateFieldWork);

location = {'SWEET','WELL3','WELL4','WELL4A','RISE'};

for iLoc = 1:length(location)

    clear n2o_ugN_L_dATM 
    ii = 0;
    for iFieldDate = 1:numel(dateFieldWork)
        if dateFieldWork(iFieldDate) > rise.dte(1)
            ind_triplicate = find(T_n2o.sampleDate == dateFieldWork(iFieldDate)...
                       & ~isnan(T_n2o.ind_triplicate) ... skip air sample
                       & strcmpi(T_n2o.location,location{iLoc}));
        
            ind_atm = find(T_n2o.sampleDate == dateFieldWork(iFieldDate)...
                       & isnan(T_n2o.ind_triplicate));
        
            ii = ii + 1;
    
            
            if isempty(ind_atm)
                ind_atm_all = find(isnan(T_n2o.ind_triplicate));
                [~, ind_atm_closest] = min(abs...
                    (T_n2o.sampleDate(ind_atm_all) - dateFieldWork(iFieldDate)));
                
                fprintf('no ATM on %s, will take the closet ATM N2O on %s at %s\n',...
                    dateFieldWork(iFieldDate),...
                    T_n2o.sampleDate(ind_atm_all(ind_atm_closest)),...
                    T_n2o.location{ind_atm_all(ind_atm_closest)})
                n2o_ugN_L_dATM(ii,1) = mean(T_n2o.n2o_ugN_L(ind_triplicate)) ...
                              % - mean(T_sample.n2o_ugN_L(ind_atm_all(ind_atm_closest)));       
            else
                n2o_ugN_L_dATM(ii,1) = mean(T_n2o.n2o_ugN_L(ind_triplicate)) ...
                              % - mean(T_sample.n2o_ugN_L(ind_atm));       
            end
        end
    end
    eval(sprintf('T_n2o_ugN_L_dATM.%s = n2o_ugN_L_dATM',location{iLoc}))
end


%% sw-gw mixing frequency
clc
clf; fntsize = 14;

xl = [rise.dte(1) rise.dte(end)+days(1)];
% xl = [rise.dte(1) rise.dte(end)];

[commonDate,iSink,iRise] = intersect(sink.dte,rise.dte);

clear ind
ii = 0;
for i=1:length(dateFieldWork)
    if dateFieldWork(i) > sink.dte(1)
        ii = ii +1;
        [~,ind(ii)]=min(abs(dateFieldWork(i)-commonDate ));
    end
end

nrow = 6;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'.k','linewidth',2); hold on
% plot(rise.dte,rise.data,'.b','linewidth',2);
% scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

% xlabel('Time');
ylabel('Flowrate (ft^3/s)');

% legend("Sink","Rise","Sampling",'location','northwest')
% legend box off
set(gca,'fontsize',fntsize,'XTickLabel','')

xlim(xl)

xl = xlim;

% yyaxis right
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-r','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


% grid on
% ------------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on

ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
ylim([-500 500])

%%
clc
clf; fntsize = 14;

xl = [rise.dte(1) rise.dte(end)+days(1)];
% xl = [rise.dte(1) rise.dte(end)];

[commonDate,iSink,iRise] = intersect(sink.dte,rise.dte);

clear ind
ii = 0;
for i=1:length(dateFieldWork)
    if dateFieldWork(i) > sink.dte(1)
        ii = ii +1;
        [~,ind(ii)]=min(abs(dateFieldWork(i)-commonDate ));
    end
end

nrow = 18;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'--k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);
scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

% xlabel('Time');
ylabel('Flowrate (ft^3/s)');

% legend("Sink","Rise","Sampling",'location','northwest')
% legend box off
set(gca,'fontsize',fntsize,'XTickLabel','')

xlim(xl)
xl = xlim;

% yyaxis right
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-r','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


% grid on
% ------------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('\Delta Q (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-12),T_n2o_ugN_L_dATM.SWEET,barwdith ,'facealpha',0.5,'facecolor','k','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-6),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 

ylabel('N_2O');% [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')

% legend("Sink - Rise","sampling","",location{2:4},'location','northwest')
% legend box off

set(gca,'fontsize',fntsize)

% ------------------------------------
% h3 = subplot(nrow ,1,7:9)
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% % scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% ax = gca;
% % xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% % yp = ax.YLim([1 1 2 2]);  
% % patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% fctr = 0.85;
% % text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% % text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)
% 
% 
% yyaxis right
% ii = 0;
% 
% % bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% % ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters
% 
% barwdith = 0.1;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-6),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% % yl = ylim;
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([0.15 0.24])
% 
% ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters
% 
% ax = gca;
% 
% ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
% % xlabel('Time');
% set(gcf,'color','w')
% 
% 
% set(gca,'fontsize',fntsize)


% ------------------------------------
subplot(nrow,1,7:9)
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('\Delta Q');% (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)



yyaxis right

barwdith = 0.1;

bar(T_DO.dateFieldWork+hours(-6),T_DO.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_DO.dateFieldWork+hours(0),T_DO.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
bar(T_DO.dateFieldWork+hours(6),T_DO.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 

bar(T_DO.dateFieldWork+hours(12),T_DO.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 

ylabel('DO');% [mg/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';

ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')
set(gca,'fontsize',fntsize,'XTickLabel','')

set(gca,'fontsize',fntsize)

% ------------------------------------
subplot(nrow,1,10:12)
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('\Delta Q');% (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)



yyaxis right

barwdith = 0.1;

bar(T_DO.dateFieldWork+hours(-6),T_DOC.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_DO.dateFieldWork+hours(0),T_DOC.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
bar(T_DO.dateFieldWork+hours(6),T_DOC.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
bar(T_DO.dateFieldWork+hours(12),T_DOC.RISE/30,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 

ylabel('DOC');% [mg C/L] '); % 1 cfs = 28.3168 liters

ax = gca;
ax.XTick = xl(1):days(5):xl(2)
ax.YColor = 'k';

% xlabel('Time');
set(gcf,'color','w')
set(gca,'fontsize',fntsize,'XTickLabel','')

set(gca,'fontsize',fntsize)

% ------------------------------------
% subplot(nrow,1,13:15)
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% ax = gca;
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% 
% 
% 
% yyaxis right
% 
% barwdith = 0.1;
% bar(T_ORP.dateFieldWork+hours(-6),T_ORP.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_ORP.dateFieldWork+hours(0),T_ORP.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_ORP.dateFieldWork+hours(6),T_ORP.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% 
% ylabel('ORP [mv] '); % 1 cfs = 28.3168 liters
% 
% ax = gca;
% 
% ax.YColor = 'k';
% 
% % xlabel('Time');
% set(gcf,'color','w')
% set(gca,'fontsize',fntsize,'XTickLabel','')
% 
% set(gca,'fontsize',fntsize)


% ------------------------------------
subplot(nrow,1,13:15)
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('\Delta Q');% (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)


% axtop = gca;


yyaxis right

barwdith = 0.1;
bar(T_no3.dateFieldWork+hours(-6),T_no3.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(0),T_no3.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(6),T_no3.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(12),T_no3.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 

ylabel('NO3');% [ugN/L] '); % 1 cfs = 28.3168 liters

ax = gca;

% ax.XTick = union([xl(1):days(5):xl(2)],dateFieldWork-timeofday(dateFieldWork))
ax.XTick = dateFieldWork-timeofday(dateFieldWork);
xtickangle(45);
ax.YColor = 'k';

% xlabel('Time');
set(gcf,'color','w')
% set(gca,'fontsize',fntsize,'XTickLabel','')

set(gca,'fontsize',fntsize)

% ------------------------------------
subplot(nrow,1,16:18)
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('\Delta Q');% (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)


% axtop = gca;


yyaxis right

barwdith = 0.1;
% table2array(T_fdom_b)
T_data = T_nh4;
% T_data(:,2:end) = (T_fdom_b(:,2:end) + T_fdom_t(:,2:end) + T_fdom_a(:,2:end));% + T_fdom_c(:,2:end));%./T_fdom(:,2:end);%.*T_DOC(:,2:end);
bar(T_no3.dateFieldWork+hours(-6),T_data.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(0),T_data.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(6),T_data.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
bar(T_no3.dateFieldWork+hours(12),T_data.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 

ylabel('NH_4');% [ugN/L] '); % 1 cfs = 28.3168 liters

ax = gca;

% ax.XTick = union([xl(1):days(5):xl(2)],dateFieldWork-timeofday(dateFieldWork))
ax.XTick = dateFieldWork-timeofday(dateFieldWork);
xtickangle(45);
ax.YColor = 'k';

% xlabel('Time');
set(gcf,'color','w')
% set(gca,'fontsize',fntsize,'XTickLabel','')

set(gca,'fontsize',fntsize)

% ------------------------------------

% subplot(nrow,1,19:21)
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% ax = gca;
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% 
% 
% 
% yyaxis right
% 
% barwdith = 0.1;
% bar(T_SpC.dateFieldWork+hours(-6),T_SpC.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_SpC.dateFieldWork+hours(0),T_SpC.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_SpC.dateFieldWork+hours(6),T_SpC.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% 
% ylabel('SpC [uS/cm] '); % 1 cfs = 28.3168 liters
% =============================

% ------------------------------------
% h2 = subplot(nrow ,1,16:18);
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% % scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% ax = gca;
% % xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% % yp = ax.YLim([1 1 2 2]);  
% % patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% fctr = 0.85;
% % text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% % text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)
% 
% 
% yyaxis right
% ii = 0;
% 
% % bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% % ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters
% 
% barwdith = 0.1;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-6),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% % yl = ylim;
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% % ylim([-0.4 2])
% % ylim([-0.5 0.8])
% ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters
% 
% ax = gca;
% 
% ax.YColor = 'k';
% ax.XTick = xl(1):days(5):xl(2)
% set(gca,'fontsize',fntsize,'XTickLabel','')
% % xlabel('Time');
% set(gcf,'color','w')
% 
% 
% set(gca,'fontsize',fntsize)
% 
% 
% % legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
% legend("Sink - Rise","","N2O",'location','northwest')


% ------------------------------------
% h3 = subplot(nrow,1,19:21);
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% % scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% ax = gca;
% % xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% % yp = ax.YLim([1 1 2 2]);  
% % patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% fctr = 0.85;
% % text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% % text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)
% 
% 
% yyaxis right
% ii = 0;
% 
% % bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% % ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters
% 
% barwdith = 0.1;
% bar(T_16s.dateFieldWork+hours(-6),T_Nir_Nos.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_16s.dateFieldWork+hours(0),T_Nir_Nos.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% % yl = ylim;
% % bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% % ylim([-0.4 2])
% 
% ylabel('Nir/Nos'); % 1 cfs = 28.3168 liters
% 
% ax = gca;
% 
% ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
% ax.XTick = xl(1):days(5):xl(2)
% % xlabel('Time');
% set(gcf,'color','w')
% 
% 
% set(gca,'fontsize',fntsize)
% 
% 
% % legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
% legend("Sink - Rise","","Nir/Nos",'location','northwest')
% 
% 
% ax = gca;
% 
% ax.YColor = 'k';
% 
% xlabel('Time');
% set(gcf,'color','w')
% 
% 
% set(gca,'fontsize',fntsize)
% 
% 
% ax.XTickLabelRotation = 45;
% ax.XTick = xl(1):days(5):xl(2)
% 
% %
% fig = gcf;
% fig.Units = 'inches';
% fig.PaperPosition = [0 0 8 6];
% 
% % flName = 'sink_rise_n2o.png'
% flNameFig = 'sink_rise_n2o.png';
% print(fig,flNameFig,'-dpng')


%%
clf
nrow = 15;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'-k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);
scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

% xlabel('Time');
ylabel('Flowrate (ft^3/s)');

legend("Sink","Rise","Sampling",'location','northwest')
legend box off
set(gca,'fontsize',fntsize,'XTickLabel','')


xlim(xl)
xl = xlim;

ax.XTick = xl(1):days(5):xl(2)
% grid on
% ------------------------------------
% h2 = subplot(nrow ,1,4:6);
% plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
% ax = gca;
% % xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% % yp = ax.YLim([1 1 2 2]);  
% % patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
% 
% plot([commonDate(1) commonDate(end)],[0 0],'--r'); 
% 
% 
% ylabel('Sink - Rise (ft^3/s)');
% 
% yl = ylim;
% thld = max(abs(yl));
% ylim([-thld thld])
% 
% xlim(xl)
% fctr = 0.85;
% % text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% % text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)
% 
% 
% yyaxis right
% ii = 0;
% 
% % bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% % ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters
% 
% barwdith = 0.1;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-12),T_n2o_ugN_L_dATM.SWEET,barwdith ,'facealpha',0.5,'facecolor','k','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-6),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% 
% ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters
% 
% ax = gca;
% 
% ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
% % xlabel('Time');
% set(gcf,'color','w')
% 
% legend("Sink - Rise","sampling","",location{:},'location','northwest')
% legend box off
% 
% set(gca,'fontsize',fntsize)

% ------------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-6),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])
% ylim([-0.5 0.8])
ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize,'XTickLabel','')
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
legend("Sink - Rise","","N2O",'location','northwest')


% ------------------------------------
h3 = subplot(nrow ,1,7:9);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_16s.dateFieldWork+hours(-6),T_Nir_Nos.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_16s.dateFieldWork+hours(0),T_Nir_Nos.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])

ylabel('Nir/Nos'); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
legend("Sink - Rise","","Nir/Nos",'location','northwest')


% ------------------------------------
h3 = subplot(nrow ,1,10:12);
ax = gca;

barwdith = 0.1;
ind_well = find(strcmpi(microbeName_well4,microbeName_well3(1)));

bar(date_16s_rDNA_well3+hours(-6),portionRate_well3(1,:),barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on

ylabel(sprintf('%s ',microbeName_well3{1})); % 1 cfs = 28.3168 liters

% yyaxis right

bar(date_16s_rDNA_well4+hours(0),portionRate_well4(ind_well,:),barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 


xlim(xl)
% ylabel(sprintf('%s \n (well4)',microbeName_well4{ind_well})); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


legend(location{[2 3 ]},'location','northwest')


% ------------------------------------
h3 = subplot(nrow ,1,13:15);
ax = gca;

barwdith = 0.1;
ind_well = find(strcmpi(microbeName_well3,microbeName_well4(1)));

bar(date_16s_rDNA_well3+hours(-6),portionRate_well3(ind_well,:),barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on

ylabel(sprintf('%s ',microbeName_well3{ind_well})); % 1 cfs = 28.3168 liters

% yyaxis right

bar(date_16s_rDNA_well4+hours(0),portionRate_well4(1,:),barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 


xlim(xl)
% ylabel(sprintf('%s \n(well4)',microbeName_well4{1})); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2);
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


legend(location{[2 3 ]},'location','northwest')

%% NSF (NO3-, N2O, nir/nos over time at Wells 3 & 4)
clf
nrow = 12;

% --------------------- Hydrograph ------------------------------
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'-k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);

scatter(commonDate(ind),zeros(size(rise.data(iSink(ind)))),200,'r^','linewidth',2)

yl = ylim;
plot([commonDate(ind) commonDate(ind)]',[yl(1)*ones(size(rise.data(iSink(ind)))) yl(2)*ones(size(rise.data(iSink(ind))))]',...
    'r:')
% scatter(commonDate(ind),yl(2)*ones(size(rise.data(iSink(ind)))),200,'r^','linewidth',2)
% scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)
ylim(yl)


ax = gca;
ylabel('Flowrate (ft^3/s)');

% legend("Sink (Q_S)","Rise (Q_R)","Sampling",'location','northeast')
% legend box off
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')


xl = [datetime('1-15-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local') ... 
      datetime('2-27-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local')];
xlim(xl)
xl = xlim;

ax.XTick = xl(1):days(5):xl(2);
% grid on

% ------------------------------------ NO3- ------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--k'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;


yyaxis right
ii = 0;


barwdith = 0.1;
bar(T_no3.dateFieldWork+hours(-4),T_no3.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on;
bar(T_no3.dateFieldWork+hours(4),T_no3.WELL4,barwdith ,'facealpha',0.5,'facecolor','M','edgealpha',1); hold on;

ylabel('NO_3^- [ug N/L] '); % 1 cfs = 28.3168 liters

% legend("Sink - Rise (Q_S - Q_R)","","WELL3","WELL4",'location','northeast')
% legend box off

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')


% ------------------------------ N2O ------------------------------------
h3 = subplot(nrow ,1,7:9);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on

ax = gca;
    
plot([commonDate(1) commonDate(end)],[0 0],'--k'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;


yyaxis right
ii = 0;


barwdith = 0.1;


bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-4),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(4),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','m','edgealpha',1); 
plot([commonDate(1) commonDate(end)],mean(T_n2o.n2o_ugN_L(ind_atm_all(ind_atm_closest)))*ones(1,2),'-r'); 

ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')
% xlabel('Time');
% set(gcf,'color','w')


% set(gca,'fontsize',fntsize)


% legend("Sink - Rise (Q_S - Q_R)","","WELL3","WELL4","Atmospheric N_2O",'location','northeast')
% legend box off


% ----------------------------- Nir/Nos -----------------------------------
h4 = subplot(nrow ,1,10:12);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--k'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_16s.dateFieldWork+hours(-6),T_Nir_Nos.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_16s.dateFieldWork+hours(0),T_Nir_Nos.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
ylim([0 20])

ylabel('Nir/Nos'); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","","Nir/Nos (WELL3)",'location','northeast')
% legend box off


% ----------------------------- 

% set(gca,'XTickLabel','')
ax.XTickLabelRotation = 30;


%
fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 8 12];


flNameFig = 'sink_rise_n2o.png';
print(fig,flNameFig,'-dpng')

'done'

%% NCKRI (N2O over time at Wells 3 & 4
clf
nrow = 9;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'-k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);
scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
ylabel('Flowrate (ft^3/s)');

legend("Sink (Q_S)","Rise (Q_R)","Sampling",'location','northeast')
legend box off
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')


xl = [datetime('1-15-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local') ... 
      datetime('2-27-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local')];
xlim(xl)
xl = xlim;

ax.XTick = xl(1):days(5):xl(2);
% grid on

% ------------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_no3.dateFieldWork+hours(-4),T_no3.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on;
bar(T_no3.dateFieldWork+hours(4),T_no3.WELL4,barwdith ,'facealpha',0.5,'facecolor','M','edgealpha',1); hold on;
bar(T_no3.dateFieldWork+hours(12),T_no3.RISE,barwdith ,'facealpha',0.5,'facecolor','G','edgealpha',1); hold on;

% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(3),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','m','edgealpha',1); 

% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])
% ylim([-0.5 0.8])
ylabel('NO_3^- [ug N/L] '); % 1 cfs = 28.3168 liters


ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')
set(gcf,'color','w')


% set(gca,'fontsize',fntsize)


% legend("Sink - Rise (Q_S - Q_R)","",location{[2 3 4]},'location','northeast')
% legend box off
% legend("NO3-","N2O",'location','northwest')



% ------------------------------------
h3 = subplot(nrow ,1,7:9);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
% bar(T_no3.dateFieldWork+hours(-3),T_no3.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on;
% bar(T_no3.dateFieldWork+hours(3),T_no3.RISE,barwdith ,'facealpha',0.5,'facecolor','M','edgealpha',1); hold on;






bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-4),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
plot([commonDate(1) commonDate(end)],mean(T_n2o.n2o_ugN_L(ind_atm_all(ind_atm_closest)))*ones(1,2),'-r'); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(4),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','m','edgealpha',1); 

bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','G','edgealpha',1); 

% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])
% ylim([-0.5 0.8])
ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize)
% set(gca,'XTickLabel','')
xlabel('Time');
set(gcf,'color','w')


% set(gca,'fontsize',fntsize)


legend("Sink - Rise (Q_S - Q_R)","","WELL3","Atmospheric N_2O","WELL4",'location','northeast')
legend box off
% legend("NO3-","N2O",'location','northwest')


% set(gca,'fontsize',fntsize)


ax.XTickLabelRotation = 30;


%
fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 9];

% flName = 'sink_rise_n2o.png'
flNameFig = 'sink_rise_n2o.png';
print(fig,flNameFig,'-dpng')
%% NCKRI (Water Isotope data)
clf; clc
scatter(T_d18O.WELL3,T_dD.WELL3,140,'rs','filled'); hold on
text(T_d18O.WELL3+.02,T_dD.WELL3,num2str([1:length(T_d18O.WELL3)]')); hold on
scatter(T_d18O.WELL4,T_dD.WELL4,150,'g^','filled')
text(T_d18O.WELL4+0.02,T_dD.WELL4,num2str([1:length(T_d18O.WELL3)]')); hold on
scatter(T_d18O.RISE,T_dD.RISE,100,'ko','filled')
text(T_d18O.RISE+.02,T_dD.RISE,num2str([1:length(T_d18O.WELL3)]')); hold on

grid on
[h,icon] = legend('Well 3','Well 4','Rise','location','northwest','fontsize',20);

% legend box off
% legend('location','northwest','fontsize',20)

icon(4).Children.MarkerSize = 20
icon(5).Children.MarkerSize = 20
icon(6).Children.MarkerSize = 20


xlabel('\Delta^{18}O [‰ VSMOW]')
ylabel('\Delta^{2}H [‰ VSMOW]')
ax =gca;
ax.FontSize = 20;

fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 6];


flNameFig = 'isotope.png';
print(fig,flNameFig,'-dpng')

%% hydrograph + N2O + Nir/Nos
clf

fntsize = 10;
nrow = 9;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'-k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);
scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
ylabel('Flowrate [ft^3/s]');

legend("Sink (Q_S)","Rise (Q_R)","Sampling",'location','northeast')
legend box off
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')


xl = [datetime('1-15-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local') ... 
      datetime('2-27-2025 00:00:00','InputFormat','MM-dd-yyyy HH:mm:ss','TimeZone','local')];
xlim(xl)
xl = xlim;

ax.XTick = xl(1):days(5):xl(2);
% grid on



% ------------------------------------
h3 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--k'); 


ylabel('\Delta Q [ft^3/s]');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
% bar(T_no3.dateFieldWork+hours(-3),T_no3.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); hold on;
% bar(T_no3.dateFieldWork+hours(3),T_no3.RISE,barwdith ,'facealpha',0.5,'facecolor','M','edgealpha',1); hold on;






bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-4),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(4),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','m','edgealpha',1); 
plot([commonDate(1) commonDate(end)],mean(T_n2o.n2o_ugN_L(ind_atm_all(ind_atm_closest)))*ones(1,2),'-r'); 

% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','G','edgealpha',1); 

% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(0),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])
% ylim([-0.5 0.8])
ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize)
set(gca,'XTickLabel','')
% 
set(gcf,'color','w')


% set(gca,'fontsize',fntsize)


legend("\Delta Q (=Q_S - Q_R)","","N2O (Well 3)","N2O (Well 4)","Atmospheric N2O",'location','northeast')
legend box off
% legend("NO3-","N2O",'location','northwest')


% set(gca,'fontsize',fntsize)





% ------------------------------------
h3 = subplot(nrow ,1,7:9);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--k'); 


ylabel('\Delta Q [ft^3/s]');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_16s.dateFieldWork+hours(-6),T_Nir_Nos.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
% bar(T_16s.dateFieldWork+hours(0),T_Nir_Nos.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])

xlabel('Time');
ylabel('Nir/Nos'); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
legend("\Delta Q (=Q_S - Q_R)","","Nir/Nos (WELL3)",'location','northeast')
legend box off

ax.XTickLabelRotation = 30;


%
fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 6];

% flName = 'sink_rise_n2o.png'
flNameFig = 'sink_rise_n2o.png';
print(fig,flNameFig,'-dpng')

%% Cham's seminar
clf
nrow = 9;
h1 = subplot(nrow ,1,1:3);

plot(sink.dte,sink.data,'-k','linewidth',2); hold on
plot(rise.dte,rise.data,'--b','linewidth',2);
scatter(commonDate(ind),rise.data(iSink(ind)),200,'r^','linewidth',2)

ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)

% xlabel('Time');
ylabel('Flowrate (ft^3/s)');

legend("Sink","Rise","Sampling",'location','northwest')
legend box off
set(gca,'fontsize',fntsize,'XTickLabel','')


xlim(xl)
xl = xlim;

ax.XTick = xl(1):days(5):xl(2)


% ------------------------------------
h2 = subplot(nrow ,1,4:6);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(-4),T_n2o_ugN_L_dATM.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(4),T_n2o_ugN_L_dATM.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])
% ylim([-0.5 0.8])
ylabel('N_2O [ug N/L] '); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
ax.XTick = xl(1):days(5):xl(2)
set(gca,'fontsize',fntsize,'XTickLabel','')
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
legend("Sink - Rise","","Well 3","Well 4",'location','northeast')
legend box off


% ------------------------------------
h3 = subplot(nrow ,1,7:9);
plot(commonDate,sink.data(iSink)-rise.data(iRise),'-k','linewidth',2); hold on
% scatter(commonDate(ind),sink.data(iSink(ind))-rise.data(iSink(ind)),200,'r^','linewidth',2)
ax = gca;
% xp = [t_tracingStart t_tracingEnd t_tracingEnd t_tracingStart];
% yp = ax.YLim([1 1 2 2]);  
% patch(xp,yp,[.5 .5 .5],'FaceAlpha',0.1)
    
plot([commonDate(1) commonDate(end)],[0 0],'--r'); 


ylabel('Sink - Rise (ft^3/s)');

yl = ylim;
thld = max(abs(yl));
ylim([-thld thld])

xlim(xl)
fctr = 0.85;
% text(min(xl)+days(1),thld*fctr,'recharging','fontsize',fntsize)
% text(min(xl)+days(1),-thld*fctr,'gaining','fontsize',fntsize)


yyaxis right
ii = 0;

% bar(commonDate(ind),n2o_rise(:).*rise.data(iRise(ind))*28.3168/1000,1,'facealpha',0.5,'facecolor','b'); 
% ylabel('N_2O [ug N]'); % 1 cfs = 28.3168 liters

barwdith = 0.1;
bar(T_16s.dateFieldWork+hours(-4),T_Nir_Nos.WELL3,barwdith ,'facealpha',0.5,'facecolor','c','edgealpha',1); 
bar(T_16s.dateFieldWork+hours(4),T_Nir_Nos.WELL4,barwdith ,'facealpha',0.5,'facecolor','r','edgealpha',1); 
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(6),T_n2o_ugN_L_dATM.WELL4A,barwdith ,'facealpha',0.5,'facecolor','y','edgealpha',1); 
% yl = ylim;
% bar(T_n2o_ugN_L_dATM.dateFieldWork+hours(12),T_n2o_ugN_L_dATM.RISE,barwdith ,'facealpha',0.5,'facecolor','b','edgealpha',1); 
% ylim([-0.4 2])

ylabel('Nir/Nos'); % 1 cfs = 28.3168 liters

ax = gca;

ax.YColor = 'k';
% set(gca,'fontsize',fntsize,'XTickLabel','')
ax.XTick = xl(1):days(5):xl(2)
% xlabel('Time');
set(gcf,'color','w')


set(gca,'fontsize',fntsize)


% legend("Sink - Rise","",location{[2 3 4]},'location','northwest')
% legend("Sink - Rise","","Nir/Nos",'location','northwest')
legend("Sink - Rise","","Well 3","Well 4",'location','northeast')
legend box off

set(gca,'fontsize',fntsize)


% legend(location{[2 3 ]},'location','northwest')
ax.XTickLabelRotation = 30;


%
fig = gcf;
fig.Units = 'inches';
fig.PaperPosition = [0 0 7 9];

% flName = 'sink_rise_n2o.png'
flNameFig = 'sink_rise_n2o.png';
print(fig,flNameFig,'-dpng')

