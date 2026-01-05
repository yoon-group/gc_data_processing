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
% Naming rule of the sample data files should be as following:
% 
%  (text 1)_(text 2)_(text 3)_(text 4)_(text 5).CSV
% 
%  where: 
%       text 1: location (e.g. RISE)
%       text 2: sampling date (e.g. 20241004)
%       text 3: tripcliate index (e.g. 3) OR ATM (e.g. ATM)
%       text 4: dilution factor (e.g. 2DF). If not diluted, skip.
%       text 5: temperature (e.g. 24.9T)
% 
% An example:
% 
%       RISE_20241004_60ML_8MIN_3_2DF_24.9T.CSV
% 
% A user should provide the value of volumes of aqueous phase (Vaq) and gas 
% phase (Vg) in the HS protocol. If not provided this code assumes as:
% 
    % Vaq : Vg = 9 : 1 
% 
% (Madison Flint) For the protocol study, we can also read Vaq and Vg from 
% file names. 
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
% 
% A user should specify the expected arrival time to Electro Capture 
% Detector (ECD) according to gas type to be detected. For N2O, our GC
% works with:
%           
    expcted_arrival_time = 1.8; % [min] 
%
% Once run, the saved result file will be used to identify new data
% which have not been analyzed, and do the additional anaysis only for the
% new data. 
% =========================================================================
%% spurious samples


%% directory setting
dir_.archive = sprintf('/home/public/gcData/');

dir_.results = sprintf('%s/processed_hsExtract/',pwd);
if ~exist(dir_.results,'dir'); mkdir(dir_.results); end

dir_.fig = [dir_.results 'fig/sample/']; 
if ~exist(dir_.fig,'dir'); mkdir(dir_.fig); end

dir_.data = sprintf('%shsExtract/',dir_.archive);

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
if exist(flNameOut,'file')
    load(flNameOut,'T_sample');
    flNameList = setdiff(flNameList,T_sample.flNameList);
end
nSample = length(flNameList);

%% location, sample date, temperature, dilution factor, triplicate index
location = cell(nSample,1);
ind_triplicate = nan(nSample,1);
dilutionFactor = ones(nSample,1);
tempK = nan(nSample,1);
bottleSize = nan(nSample,1);
shakeTime = nan(nSample,1);
suspicious = false(nSample,1);
for iSample = 1:nSample

    note = textscan(flNameList{iSample}(1:end-4), '%s','delimiter','_' ); 
    note = note{1};

    location{iSample} = note{1};
    sampleDate(iSample,1) = datetime(note{2},'InputFormat','yyyyMMdd','TimeZone','local');
    
    if numel(note)>2 
        for iText = 3:numel(note)
            if ~contains(upper(note{iText}),'ATM')
                if ~isempty(str2num(note{iText}))
                    ind_triplicate(iSample) = str2num(note{iText});
                end
                
            end

            if contains(upper(note{iText}),'SUSPICIOUS')
                suspicious(iSample) = true;
            end
            if contains(upper(note{iText}),'DF')
                dilutionFactor(iSample) = str2num(erase(note{iText},'DF'));
            end
            if contains(upper(note{iText}),'T')
                tempK(iSample) = str2num(erase(note{iText},'T')) + 273.15;
            end
            if contains(upper(note{iText}),'MIN')
                shakeTime(iSample) = str2num(erase(note{iText},'MIN'));
            end
            if contains(upper(note{iText}),'ML')
                tmp = str2num(erase(note{iText},'ML'));
                switch tmp
                    case 150
                        bottleSize(iSample) = 157;
                    case 60
                        bottleSize(iSample) = 70;
                    case 20
                        bottleSize(iSample) = 24;
                    case 7.5
                        bottleSize(iSample) = 7.84;
                end
            end
        end
    end
end

%% peak measuring
peakH = nan(nSample,1);
for iSample = 1:length(flNameList)
    fprintf('(%d/%d) %s --> processing -->',iSample,length(flNameList),flNameList{iSample}(1:end-4))
    peakH(iSample) = peakMeasure(flNameList{iSample},expcted_arrival_time,dir_);    
    fprintf('  completed \n')
end

%% load calibration data 
flNameSTD = sprintf('%s/processed/STD.mat',pwd);
load(flNameSTD,'T_STD','slope','offset')

%% N2O concentration estimation
bottleVol = bottleSize/1000; %[L]
Vaq = bottleVol*9/10; %[L] 
Vg =  bottleVol*1/10; %[L] 
valveVol = 3/1000; %[L]
n2o_ppm_atm = 340/1000; %[ppm]

% Vaq = 0.54; Vg = 0.06; % [L] depending on your HS protocol

if nSample 
    logK0 = A + (B./tempK) + C*(log(tempK));
    K0 = exp(logK0); % Henry's law constant [mol/L/atm]
    
    n2o_ppm_measured = (peakH*slope + offset).*dilutionFactor;
    
    n2o_ppm = (n2o_ppm_measured .* (Vg+valveVol) - (valveVol*n2o_ppm_atm))./Vg;


    n2o_atm = n2o_ppm * 1e-6; % Fractional concentration of N2O in the gas 
                              % phase. By mulitiplying the total pressure 
                              % of 1 atmospheres (atm), the value also 
                              % expresses a partial pressure [atm].
                              % (Henry's law: naq/Vaq = KO * n2o_atm)
    
    ng = (n2o_atm .* Vg)./(R*tempK); % N2O moles in gas phase (from ideal gas law)
    naq = K0.*n2o_atm.*Vaq; % dissolved N2O moles in aqueous phase (from Henry's law)
    
    n2o_mol_L = (ng+naq)./Vaq; % Total N2O concentration in mol/L
    n2o_nmol_L = n2o_mol_L * 1e+9;
    n2o_umol_L = n2o_mol_L * 1e+6;
    
    n2o_ng_L = n2o_nmol_L * amuN2O;
    n2o_ug_L = n2o_umol_L * amuN2O;
    
    n2o_molN_L = n2o_mol_L * 2;
    n2o_gN_L = n2o_molN_L * amuN;
    n2o_mgN_L = n2o_gN_L * 1e+3;
    n2o_ugN_L = n2o_mgN_L * 1e+3;
end

%% add newly processed data & save
if nSample
    T_new = table(flNameList,location,suspicious,bottleSize,shakeTime,ind_triplicate,...
        dilutionFactor,tempK,peakH,n2o_ppm,n2o_molN_L,n2o_ugN_L,sampleDate);
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
