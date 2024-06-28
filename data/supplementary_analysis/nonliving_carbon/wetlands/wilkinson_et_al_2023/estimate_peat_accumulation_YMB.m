%==========================================================================
% Script to run Monte Carlo simulation for examining the carbon sink
%   function of northern peatlands given: (i) variability of post-fire NEE;
%   (ii) "recovered" peatland NEE; (iii) time for recovery; (iv) and fire
%   return interval.
%
% Script includes figures generated from the Monte Carlo simulations used 
%   in the manuscript submitted to Nature Climate Change:
%   NCLIM-22071425B
%==========================================================================

%% Examine NEE data

% Using current directory, so assumes input files are in the same folder as
%   the 'Wilkinson_et_al_NCC.m' file.
folder = cd;
%--> Excel spreasheet which contains annual NEE estimates from a literature
%       search of northern non-permafrost peatlands.
file = 'Annual_NEE.xlsx';

% Read in burned peatland site NEE
data = readtable(fullfile(folder,file),'Sheet','burned');

% These are annual estimates of NEE including winter fluxes
NEE.burned = data.Estimated_CO2_Annual_g_C_flux;

% Read in northern non-permafrost peatland site NEE
data = readtable(fullfile(folder,file),'Sheet','not_burned');

% These are annual estimates of NEE including winter fluxes (measured or
%   modelled).
NEE.raw = data.Estimated_CO2_Annual_g_C_flux;

% Generate indexes for different peatland states
idx.degradation = ~strcmpi(data.Drainage_type,'None') & ~strcmpi(data.Drainage_type,'intact');
idx.drained = strcmpi(data.Drainage_type,'drained') | strcmpi(data.Drainage_type,'degraded'); % Drained for agriculture, forestry, or horticulture
idx.rewetted = strcmpi(data.Drainage_type,'rewetted') | strcmpi(data.Drainage_type,'restored'); % Rewetted with or without active veg restoration (e.g. MLTT)
idx.notNAN = ~isnan(NEE.raw);

NEE.recovered = NEE.raw(~idx.degradation & idx.notNAN);
NEE.drained = NEE.raw(idx.drained & idx.notNAN);
NEE.rewetted = NEE.raw(idx.rewetted & idx.notNAN);

% Used 'isoutlier to objectively remove outliers from data sets using a
%   non-parametric method.
NEE.recovered = NEE.recovered(~isoutlier(NEE.recovered,'median'));
NEE.burned = NEE.burned(~isoutlier(NEE.burned,'median'));
NEE.drained = NEE.drained(~isoutlier(NEE.drained,'median'));
NEE.rewetted = NEE.rewetted(~isoutlier(NEE.rewetted,'median'));

% Array containing annual NEE data from the compiled literature data set. 
X = [NEE.burned; NEE.recovered; NEE.drained; NEE.rewetted];
grp = [ones(length(NEE.burned),1); 2.*ones(length(NEE.recovered),1); 3.*ones(length(NEE.drained),1); 4.*ones(length(NEE.rewetted),1)];

% Visualize difference in NEE between burned, recovered, drained and
%   rewetted from the compiled literature data set (figure not part of the
%   manuscript).
% boxplot(X,grp)

% Distributions
%--> Original was normal distribution for all NEE groups
%--> [mean(NEE.recovered) std(NEE.recovered)./sqrt(length(NEE.recovered))]
%--> [mean(NEE.burned) std(NEE.burned)./sqrt(length(NEE.burned))]
%--> [mean(NEE.drained) std(NEE.drained)./sqrt(length(NEE.drained))]
%--> [mean(NEE.rewetted) std(NEE.rewetted)./sqrt(length(NEE.rewetted))]
%--> recovered: EpsilonSkewNormal, theta=-40.91; sigma=56.42; epsilon=0.261
%       but 'icdf' does not support EpsilonSkewNormal
%               Normal, mu=-62.6; sigma=57.8; SE=4.3
%--> drained: Normal, mu=191.8; sigma=249.5 SE=27.22
%--> rewetted: Normal, mu=-5.72; sigma=84.7; SE=12.88
%--> burned: Normal, mu=71.43; sigma=53.6; SE=14.87


%% Set values common to each MC simulation

% Number of parameter sets to generate for the Monte Carlo simulation
NU = 10000;

%% Annual g C-CH4 m^-2 yr^-1 from Abdalla et al., 2016
%--> DOI: https://doi.org/10.1002/ece3.2469

ch4.type = {'natural','drained','restored'};

ch4.drained = [2.1 0 0.7 -0.1 0.9 3.3 -0.2 0.2 0.2 0.1 0 0 0 0 0 0 0 0 -0.1 -0.1 0 1 0.7 -0.1 0.2 -0.1 -0.1 -0.2 -0.1 5.1 1.4 0.4 0 0 1.1 1 1.1 1 1 1.2 7.4 2.3 0.5 0.1 7.9 2.4 1.1 -0.2 0.9 2.4 0.2 1.2 3.3 1.8 1.3 0.2 1.7 4.0 0.7 0.7 0 0.2 0.8]';
ch4.restored = [0.1 0.9 0 1.7 0.7 0.6 9 31.8 0.9 1.6 3.5 6.9 2.7 0 0 0.4 15.5 1.4 5.9 1.2 47.8 74.7 111.4 0.2 0 0 0.1 2.2 24.6 -0.1 0.1 0.1 -0.1 17.9 8.2 9.9 5.3 14.1]';
ch4.natural = [4 0.5 6.2 31.8 4.2 2.7 0.6 0.5 0 3.8 3.2 6.7 9 3.9 0.2 0.8 2.3 9 26.7 41.1 0.3 15 154.1 102.7 1.5 11.8 0.9 2.4 1.9 8.4 4 5.3 2.7 1 0.6 5.1 4.1 18.1 16.3 11 0.1 20.3 18.3 93.3 30 29.2 7.9 4.7 4.9 3 8 0.9 0.2 0 2.8 4.3 1.5 1.8 1.1 4.4 3.7 1.9 6.2 12.4 11.5 26 6.9 16.4 24.7 1 1.9 19.4 4.1 2.9 4.9 15 9.4 1.3 0.3 0.3 0.5 3.5 3.0 0.1 0.1 5 0 3.7 53.7 18.8 1.8 2.8 2.2 6.6 17.7 20.2 22.6 31 6.3 2.8 51.8 10.4 8.9 10.1 4 2.4 7.1 5.8]';

% Re-sample to get empirical distribution of same size as Monte Carlo data
for i=1:3 % Loop over ch4.type
    idx.resample = ceil(length(ch4.(ch4.type{i})).*rand(NU,1));
    ch4.resample.(ch4.type{i}) = ch4.(ch4.type{i})(idx.resample);
end


%% Monte Carlo simulation of NEE for undisturbed northern non-permafrost peatlands

% Using current directory, so assumes \safe_R1.1\ is a subfolder in the 
%   same folder as the 'Wilkinson_et_al_NCC.m' file.
basefolder = cd;

% 'the Software': SAFE V1.1  - Sensitivity Analysis For Everybody Matlab/Octave 
%   toolbox developed by Francesca Pianosi, Thorsten Wagener and Fanny Sarrazin 
%   and owned by the University of Bristol.
% https://safetoolbox.github.io/
% https://doi.org/10.1016/j.envsoft.2015.04.009
addpath(fullfile(basefolder,'safe_R1.1'))
addpath(fullfile(basefolder,'safe_R1.1/PAWN/'))
addpath(fullfile(basefolder,'safe_R1.1/VBSA'))
addpath(fullfile(basefolder,'safe_R1.1/sampling'))
addpath(fullfile(basefolder,'safe_R1.1/util'))
addpath(fullfile(basefolder,'safe_R1.1/visualization'))

% Define input distribution and ranges:
%--> Number of parameters
M  = 6 ;
%--> Parameter distribution name
DistrFun  = {'norm','norm','logn','exp','unif','unif'} ;
% Parameter names
labelparams={ 'NEE-recovered','NEE-burned',...
    'C-loss-fire','FFI','years-burn-lag','years-2-recovery'} ; 

%--> recovered: Normal, mu=-62.6; sigma=57.8; SE=4.3
%--> burned: Normal, mu=71.43; sigma=53.6; SE=14.87

% Distribution parameters:
DistrPar=[{[-62.6,57.8]},{[71.4 53.6]},...
    {[0.586887 0.907298]},{[0.3448]},{[1 10]},{[11, 60]}];

% Sample parameter space using the resampling strategy proposed by 
% (Saltelli, 2008; for reference and more details, see help of functions
% vbsa_resampling and vbsa_indices) 

% Create parameter samples to estimate the unconditional CDF:
%--> Sample parameter space using Latin Hypercube method
Xu = AAT_sampling('lhs',M,DistrFun,DistrPar,NU);
%--> Modify FFI
Xu(:,4) = 100./Xu(:,4);

% Generate average NEE estimates using the C_balance_model and Monte Carlo
%   parameter inputs.
Yu = model_execution('C_balance_model',Xu);


%% Monte Carlo simulation of NEE using drained northern non-permafrost peatlands

% Using current directory, so assumes \safe_R1.1\ is a subfolder in the 
%   same folder as the 'Wilkinson_et_al_NCC.m' file.
basefolder = cd;

% 'the Software': SAFE V1.1  - Sensitivity Analysis For Everybody Matlab/Octave 
%   toolbox developed by Francesca Pianosi, Thorsten Wagener and Fanny Sarrazin 
%   and owned by the University of Bristol.
% ADD URL
addpath(fullfile(basefolder,'safe_R1.1'))
addpath(fullfile(basefolder,'safe_R1.1/PAWN/'))
addpath(fullfile(basefolder,'safe_R1.1/VBSA'))
addpath(fullfile(basefolder,'safe_R1.1/sampling'))
addpath(fullfile(basefolder,'safe_R1.1/util'))
addpath(fullfile(basefolder,'safe_R1.1/visualization'))

% Define input distribution and ranges:
%--> Number of parameters
M  = 6 ;
%--> Parameter distribution name
DistrFun  = {'norm','norm','logn','exp','unif','unif'} ;
% Parameter names
labelparams={ 'NEE-recovered','NEE-burned',...
    'C-loss-fire','FFI','years-burn-lag','years-2-recovery'} ; 

%--> drained: Normal, mu=191.8; sigma=249.5 SE=27.22
%--> burned: Normal, mu=71.43; sigma=53.6; SE=14.87

% Distribution parameters:
DistrPar=[{[191.8, 249.5]},{[71.4, 53.6]},...
    {[1.84561 0.845759]},{[0.3448]},{[1 10]},{[11, 60]}];

Xmin = [-300 -10 2.3 0 1 11];
Xmax = [800 200 16.8 Inf 10 60];

% Create parameter samples to estimate the unconditional CDF:
%--> Sample parameter space using Latin Hypercube method
Xu_drained = AAT_sampling('lhs',M,DistrFun,DistrPar,NU); % matrix (NU,M)
%--> Limit parameters to observed range where appropriate
for i=1:M
    idx = Xu_drained(:,i)<Xmin(i);
    Xu_drained(idx,i) = Xmin(i);
    idx = Xu_drained(:,i)>Xmax(i);
    Xu_drained(idx,i) = Xmax(i);
end

%--> Modify FFI
Xu_drained(:,4) = 100./Xu_drained(:,4);

% Generate average NEE estimates using the C_balance_model and Monte Carlo
%   parameter inputs.
Yu_drained = model_execution('C_balance_model',Xu_drained);


%% Monte Carlo simulation of NEE using restored northern non-permafrost peatlands

% Using current directory, so assumes \safe_R1.1\ is a subfolder in the 
%   same folder as the 'Wilkinson_et_al_NCC.m' file.
basefolder = cd;

% 'the Software': SAFE V1.1  - Sensitivity Analysis For Everybody Matlab/Octave 
%   toolbox developed by Francesca Pianosi, Thorsten Wagener and Fanny Sarrazin 
%   and owned by the University of Bristol.
% ADD URL
addpath(fullfile(basefolder,'safe_R1.1'))
addpath(fullfile(basefolder,'safe_R1.1/PAWN/'))
addpath(fullfile(basefolder,'safe_R1.1/VBSA'))
addpath(fullfile(basefolder,'safe_R1.1/sampling'))
addpath(fullfile(basefolder,'safe_R1.1/util'))
addpath(fullfile(basefolder,'safe_R1.1/visualization'))

% Define input distribution and ranges:
%--> Number of parameters
M  = 6 ;
%--> Parameter distribution name
DistrFun  = {'norm','norm','logn','exp','unif','unif'} ;
% Parameter names
labelparams={ 'NEE-recovered','NEE-burned',...
    'C-loss-fire','FFI','years-burn-lag','years-2-recovery'} ; 

%--> rewetted: Normal, mu=-5.72; sigma=84.7; SE=12.88
%--> burned: Normal, mu=71.43; sigma=53.6; SE=14.87

% Distribution parameters:
DistrPar=[{[-5.7,84.7]},{[71.4, 53.6]},...
    {[0.586887 0.907298]},{[0.3448]},{[1 10]},{[11, 60]}];

Xmin = [-200 -10 0.5 0 1 11];
Xmax = [200 200 10 Inf 10 60];

% Create parameter samples to estimate the unconditional CDF:
%--> Sample parameter space using Latin Hypercube method
Xu_restored = AAT_sampling('lhs',M,DistrFun,DistrPar,NU);
%--> Limit parameters to observed range where appropriate
for i=1:M
    idx = Xu_restored(:,i)<Xmin(i);
    Xu_restored(idx,i) = Xmin(i);
    idx = Xu_restored(:,i)>Xmax(i);
    Xu_restored(idx,i) = Xmax(i);
end

%--> Modify FFI
Xu_restored(:,4) = 100./Xu_restored(:,4);
% Generate average NEE estimates using the C_balance_model and Monte Carlo
%   parameter inputs.
Yu_restored = model_execution('C_balance_model',Xu_restored);

%% Generate estimate 

peat_area = 2E12; % Boreal peatland area in square meters

NCB = 0.93.*Yu;
pristine_peatlands = mean(NCB).*peat_area./(1e6)./(1e09);
pristine_peatlands 


