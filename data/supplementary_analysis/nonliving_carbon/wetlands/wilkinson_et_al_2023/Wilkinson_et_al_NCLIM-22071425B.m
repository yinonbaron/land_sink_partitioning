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


%% Figure 01 - Comparison of PDFs

figure('color','white','units','inches','position',[2 2 4 3.5])
xi = (-300:4:800)';
f = ksdensity(Yu + ch4.resample.natural, xi);
plot(xi,f,'color',[0 0 1],'linestyle','-','linewidth',2)
hold on; box on
f = ksdensity(Yu_drained + ch4.resample.drained, xi);
plot(xi,f,'color',[1 0 0],'linestyle','-','linewidth',2)
f = ksdensity(Yu_restored + ch4.resample.restored, xi);
plot(xi,f,'color',[0.5 0.5 0.5],'linestyle','-.','linewidth',2)

ylabel('Relative frequency')
xlabel('NEE + CH_4 (gC m^{-2} yr^{-1})')
set(gca,'xminortick','on','yminortick','on')


clear idx file
folder = cd;
file.recovered = 'Annual_NEE.xlsx';

% Read in recovered peatland site NEE
data = readtable(fullfile(folder,file.recovered),'Sheet','not_burned');
NEE.raw = data.Estimated_CO2_Annual_g_C_flux;

idx.degradation = ~strcmpi(data.Drainage_type,'None') & ~strcmpi(data.Drainage_type,'intact');
idx.notNAN = ~isnan(NEE.raw);

NEE.recovered = NEE.raw(~idx.degradation & idx.notNAN);
NEE.recovered = NEE.recovered(~isoutlier(NEE.recovered,'median'));
%--> Resample so that NEE and CH4 can be combined
idx.resample = ceil(length(NEE.recovered).*rand(NU,1));
NEE.resampled.natural = NEE.recovered(idx.resample);

NEE.degraded = NEE.raw(idx.degradation & idx.notNAN);
idx.resample = ceil(length(NEE.degraded).*rand(NU,1));
NEE.resampled.degraded = NEE.degraded(idx.resample);

idx.rewetted = strcmpi(data.Drainage_type,'rewetted') | strcmpi(data.Drainage_type,'resotred');
NEE.restored = NEE.raw(idx.rewetted & idx.notNAN);
idx.resample = ceil(length(NEE.restored).*rand(NU,1));
NEE.resampled.restored = NEE.restored(idx.resample);

xi = (-300:4:440)';
f = ksdensity(Xu(:,1) + ch4.resample.natural, xi);
plot(xi,f,'b:','linewidth',2)

xi = (-300:4:800)';
f = ksdensity(Xu_drained(:,1) + ch4.resample.drained, xi);
plot(xi,f,'r:','linewidth',2)

set(gca,'xlim',[-300 800])

legend('Natural','Degraded','Restored','No burn - natural','No burn - degraded')


[mean(NEE.resampled.natural+ch4.resample.natural) std(NEE.resampled.natural+ch4.resample.natural)]
[mean(Yu+ch4.resample.natural) std(Yu+ch4.resample.natural)]
[mean(Yu_drained+ch4.resample.drained) std(Yu_drained+ch4.resample.drained)]
[mean(Yu_restored+ch4.resample.restored) std(Yu_restored+ch4.resample.restored)]
[mean(Xu_drained(:,1)+ch4.resample.drained) std(Xu_drained(:,1)+ch4.resample.drained)]


%% Figure 2a - Surface plot showing influence of percent drained and FFI on median NEE+CH4

% Colour ramp adapted from IPCC report
clrSpec2 = [2 117 118; 0 143 145; 71 166 169; 116 182 185; 175 212 213;...
    237 217 181; 229 169 124; 223 139 90; 204 111 75; 174 81 76]./255;

%--> recovered: Normal, mu=-62.6; sigma=57.8; SE=4.3
%--> drained: Normal, mu=191.8; sigma=249.5 SE=27.22
%--> rewetted: Normal, mu=-5.72; sigma=84.7; SE=12.88
%--> burned: Normal, mu=71.43; sigma=53.6; SE=14.87

n_div = 50;
NCB_surf = nan(n_div,n_div);
x = nan(n_div,n_div);
y = nan(n_div,n_div);
pct_drained = linspace(0,0.49,n_div);
FFI = round(100./linspace(0,2.5,n_div));
iters = NU/2;
NEE_recovered = random('normal',-62.6, 4.3,iters,1);
NEE_drained = random('normal',191.8, 27.2,iters,1);
NEE_burned = random('normal',71.43, 14.9, iters,1);
years_burn_lag = random('unif',1,10,iters,1);
years_2_recovery = random('unif',11,60,iters,1);
C_loss_natural = random('logn',0.586887, 0.907298,iters,1);
C_loss_drained = random('logn',1.84561, 0.845759, iters,1);
for i=1:n_div % Percent drained
    for j=1:n_div % FFI
        NCB = nan(iters,1);
        
        for k=1:iters
            % Weight peat carbon loss from wildfire based on the proportion
            %   of drained peatlands.
            C_loss_fire = C_loss_natural(k).*(1-pct_drained(i)) + pct_drained(i).*C_loss_drained(k);
            % Weight the recovered NEE used in the C_balance_model based on
            %   the proportion of drained peatlands.
            NEE_for_model = NEE_recovered(k).*(1-pct_drained(i)) + pct_drained(i).*NEE_drained(k);
            % Weight the CH4 based on the proportion of drained peatlands.
            ch4.weighted = ch4.resample.natural(k).*(1-pct_drained(i)) +...
                ch4.resample.drained(k).*pct_drained(i);
            if j>1
                NCB(k,1) = C_balance_model(NEE_for_model,NEE_burned(k),C_loss_fire,...
                    FFI(j),years_burn_lag(k), years_2_recovery(k)) + ch4.weighted;
            else
                NCB(k,1) = NEE + ch4.weighted;
            end
        end
        NCB_surf(i,j) = prctile(NCB,50);
        
        % Transfer results to x-y space for surface plot
        x(i,j) = pct_drained(i);
        y(i,j) = FFI(j);
    end
end
pct_peat_area = 100./y;
peat_area = 2E12; % Boreal peatland area in square meters
NCB_GtC = NCB_surf.*peat_area./(1e6)./(1e09);


figure('color','white','units','inches','position',[2 2 6.5 3.5])
axes('position',[0.075 0.15 0.25 0.6])
xplot = nan(n_div+2,n_div+2);
xplot(2:end-1,2:end-1) = 100.*x;
yplot = nan(n_div+2,n_div+2);
yplot(2:end-1,2:end-1) = pct_peat_area;
zplot = nan(n_div+2,n_div+2);
zplot(2:end-1,2:end-1) = NCB_GtC;
zplot(zplot>0.5) = 0.5;
zplot(1,1) = -0.5;
zplot(end,end) = 0.5;
hold on

contourf(xplot,yplot,zplot,-0.5:0.1:0.5)
xlabel('Percent degraded')
ylabel('Peatland burn rate (% yr^{-1})')
set(gca,'xlim',[0 49],'ylim',[0 2.5])

colormap(clrSpec2)
[C,h] = contour3(100.*x,pct_peat_area,NCB_GtC,-0.4:0.1:0.4,'k');
clabel(C,h,'manual')
set(h,'linestyle','none')

% Add reference point of FFI=~290 (0.3448) and drained=7%
hold on; box on
plot3(7,0.3448,1,'ko','markerfacecolor','k')

text(2,2.4,0.6,'a)')

ah = gca;
ah.XMinorTick = 'on';
ah.YMinorTick = 'on';
ah.TickLength = [0.02 0.005];
ah.YAxis.TickDirection = 'both';
ah.XAxis.TickDirection = 'both';
ah.LineWidth = 1;
% set box property to off and remove background color
set(ah,'box','off')
% create new, empty axes with box but without ticks
axes('Position',get(ah,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',1);
% set original axes as active
axes(ah)


%% Figure 2b - Surface plot showing influence of FFI and smouldering.

n_div = 50;
NCB_surf = nan(n_div,n_div);
x = nan(n_div,n_div);
y = nan(n_div,n_div);
pct_drained = linspace(0.07,0.07,n_div);
C_loss_fire = linspace(0.5,10,n_div);
% FFI = round(logspace(1.6021,4,n_div));
FFI = round(100./linspace(0.01,2.5,n_div));
iters = NU/2;

NEE_recovered = random('normal',-62.6, 4.3,iters,1);
NEE_drained = random('normal',191.8, 27.2,iters,1);
NEE_burned = random('normal',71.43, 14.9, iters,1);
years_burn_lag = random('unif',1,10,iters,1);
years_2_recovery = random('unif',11,60,iters,1);
for i=1:n_div % C-loss
    for j=1:n_div % % FFI
        NCB = nan(iters,1);
        
        for k=1:iters
            NEE_for_model = NEE_recovered(k).*(1-pct_drained(i)) + pct_drained(i).*NEE_drained(k);
            ch4.weighted = ch4.resample.natural(k).*(1-pct_drained(i)) +...
                ch4.resample.drained(k).*pct_drained(i);
            NCB(k,1) = C_balance_model(NEE_for_model,NEE_burned(k),C_loss_fire(i),...
                FFI(j),years_burn_lag(k), years_2_recovery(k)) + ch4.weighted;
        end
        NCB_surf(i,j) = prctile(NCB,50);
        x(i,j) = C_loss_fire(i);
        y(i,j) = FFI(j);
    end
end
pct_peat_area = 100./y;
RCB = 120E15; % gC [equivalent to 440 Gt CO2]
peat_area = 2E12; % Boreal peatland area in square meters
NEE_noFire = 60; % gC m^-2 yr^-1
pct_RCB = 100.*(NCB_surf+NEE_noFire).*peat_area./RCB;
NCB_GtC = NCB_surf.*peat_area./(1e6)./(1e09);



axes('position',[0.4 0.15 0.25 0.6])
hold on
xplot = nan(n_div+2,n_div+2);
xplot(2:end-1,2:end-1) = x;
yplot = nan(n_div+2,n_div+2);
yplot(2:end-1,2:end-1) = pct_peat_area;
zplot = nan(n_div+2,n_div+2);
zplot(2:end-1,2:end-1) = NCB_GtC;
zplot(zplot>0.5) = 0.5;
zplot(1,1) = -0.5;
zplot(end,end) = 0.5;
% surf(xplot, yplot , zplot,'linestyle','none'); view(2)
contourf(xplot,yplot,zplot,-0.5:0.1:0.5)
xlabel('Peat carbon loss (kg C m^{-2})')
ylabel('Peatland burn rate (% yr^{-1})')
set(gca,'xlim',[0.5 10],'ylim',[0.01 2.5])
ch = colorbar('northoutside');
ch.Label.String = 'NEE + CH_4 flux (GtC yr^{-1})';
ch.Limits = [-0.2 0.5];
colormap(clrSpec2)
% colormap(cmap)
[C,h] = contour3(x,pct_peat_area,NCB_GtC,-0.4:0.1:0.5,'k');
clabel(C,h,'manual')
set(h,'linestyle','none')

hold on; box on
plot3(3.2,0.3448,1,'ko','markerfacecolor','k')

set(gca,'position',[0.4 0.15 0.25 0.6])

text(1.1,2.4,0.6,'b)')

ah = gca;
ah.XMinorTick = 'on';
ah.YMinorTick = 'on';
ah.TickLength = [0.02 0.005];
ah.YAxis.TickDirection = 'both';
ah.XAxis.TickDirection = 'both';
ah.LineWidth = 1;
% set box property to off and remove background color
set(ah,'box','off')
% create new, empty axes with box but without ticks
axes('Position',get(ah,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',1);
% set original axes as active
axes(ah)


%% Figure 2c- Surface plot showing influence of percent drained and smouldering.
%--> Using a median FFI of 290 years based on FiredPY data set for northern
%       non-permafrost peatlands.
n_div = 50;
NCB_surf = nan(n_div,n_div);
x = nan(n_div,n_div);
y = nan(n_div,n_div);
pct_drained = linspace(0,0.49,n_div);
C_loss_fire = linspace(0.5,10,n_div);
iters = NU/2;
NEE_recovered = random('normal',-62.6, 4.3,iters,1);
NEE_drained = random('normal',191.8, 27.2,iters,1);
NEE_burned = random('normal',71.43, 14.9, iters,1);
years_burn_lag = random('unif',1,10,iters,1);
years_2_recovery = random('unif',11,60,iters,1);

for i=1:n_div % Percent drained
    for j=1:n_div % C-loss
        NCB = nan(iters,1);
        
        for k=1:iters
            NEE_for_model = NEE_recovered(k).*(1-pct_drained(i)) + pct_drained(i).*NEE_drained(k);
            ch4.weighted = ch4.resample.natural(k).*(1-pct_drained(i)) +...
                ch4.resample.drained(k).*pct_drained(i);
            NCB(k,1) = C_balance_model(NEE_for_model,NEE_burned(k),C_loss_fire(j),...
                290,years_burn_lag(k), years_2_recovery(k)) + ch4.weighted;
        end
        NCB_surf(i,j) = prctile(NCB,50);
        x(i,j) = pct_drained(i);
        y(i,j) = C_loss_fire(j);
    end
end
pct_peat_area = 100./y;
RCB = 120E15; % gC [equivalent to 440 Gt CO2]
peat_area = 2E12; % Boreal peatland area in square meters
NEE_noFire = 60; % gC m^-2 yr^-1
pct_RCB = 100.*(NCB_surf+NEE_noFire).*peat_area./RCB;
NCB_GtC = NCB_surf.*peat_area./(1e6)./(1e09);


axes('position',[0.725 0.15 0.25 0.6])

hold on
xplot = nan(n_div+2,n_div+2);
xplot(2:end-1,2:end-1) = 100.*x;
yplot = nan(n_div+2,n_div+2);
yplot(2:end-1,2:end-1) = y;
zplot = nan(n_div+2,n_div+2);
zplot(2:end-1,2:end-1) = NCB_GtC;
zplot(zplot>0.5) = 0.5;
zplot(1,1) = -0.5;
zplot(end,end) = 0.5;

contourf(xplot,yplot,zplot,-0.5:0.1:0.5)
xlabel('Percent degraded')
ylabel('Peat carbon loss (kg C m^{-2})')
set(gca,'xlim',[0 49],'ylim',[0.5 10])

colormap(clrSpec2)
[C,h] = contour3(100.*x,y,NCB_GtC,-0.4:0.1:0.4,'k');
clabel(C,h,'manual')
set(h,'linestyle','none')

hold on; box on
plot3(7,3.2,1,'ko','markerfacecolor','k')

text(0.2,9.5,0.6,'c)')

ah = gca;
ah.XMinorTick = 'on';
ah.YMinorTick = 'on';
ah.TickLength = [0.02 0.005];
ah.YAxis.TickDirection = 'both';
ah.XAxis.TickDirection = 'both';
ah.LineWidth = 1;
% set box property to off and remove background color
set(ah,'box','off')
% create new, empty axes with box but without ticks
axes('Position',get(ah,'Position'),'box','on','xtick',[],'ytick',[],'linewidth',1);
% set original axes as active
axes(ah)


%% Figure 03
clear idx

peat_area = 2E12; % Boreal peatland area in square meters

% No burn
ch4.weighted = 0.93.*ch4.resample.natural + 0.07.*ch4.resample.drained;
NCB = 0.93.*NEE.resampled.natural + 0.07.*NEE.resampled.degraded + ch4.weighted;
NCB = reshape(NCB(1:30*round(NU/30)),30,round(NU/30));
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_2050 = [mean(tmp) std(tmp)];

NCB = 0.93.*NEE.resampled.natural + 0.07.*NEE.resampled.degraded + ch4.weighted;
NCB = reshape(NCB,80,NU/80);
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_2100 = [mean(tmp) std(tmp)];


% Current state (i.e. with 7% degraded)
NCB = 0.93.*Yu + 0.07.*Yu_drained + ch4.weighted;
NCB = reshape(NCB(1:30*round(NU/30)),30,round(NU/30));
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_2050(2,:) = [mean(tmp) std(tmp)];

NCB = 0.93.*Yu + 0.07.*Yu_drained + ch4.weighted;
NCB = reshape(NCB,80,NU/80);
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_2100(2,:) = [mean(tmp) std(tmp)];

% Current state with restoration
ch4.weighted = ch4.resample.natural.*(1-0.07) + ch4.resample.restored.*0.07;
NCB = 0.93.*Yu + 0.07.*Yu_restored + ch4.weighted;
NCB = reshape(NCB(1:30*round(NU/30)),30,round(NU/30));
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_2050(3,:) = [mean(tmp) std(tmp)];

NCB = 0.93.*Yu + 0.07.*Yu_restored + ch4.weighted;
NCB = reshape(NCB,80,NU/80);
tmp = mean(NCB).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_2100(3,:) = [mean(tmp) std(tmp)];


%--> With climate change scenarios
iters = 10000;
NCB_CC = nan(2000,80,2,2);
AAB = linspace(0.35,0.7,80)';
NEE_recovered = random('normal',-62.6, 4.3,5.*iters,1);
NEE_drained = random('normal',191.8, 27.2,5.*iters,1);
NEE_burned = random('normal',71.43, 14.9, 5.*iters,1);
years_burn_lag = random('unif',1,10,5.*iters,1);
years_2_recovery = random('unif',11,60,5.*iters,1);
C_loss_natural = random('logn',0.586887, 0.907298,5.*iters,1);
C_loss_drained = random('logn',1.84561, 0.845759, 5.*iters,1);
C_loss_increase = linspace(1.5,1.5,80)';
for i=1:2 % Fire severity
    for ii=1:2 % burn rate
        for j=1:80 % Years into future
            for k=1:iters
                ind = ceil(10000.*rand(1));
                C_loss_fire = C_loss_natural(ind).*0.93 + 0.07.*C_loss_drained(ind) + C_loss_increase(j)*(i-1);
                NEE.weighted = NEE_recovered(ind).*0.93 + 0.07.*NEE_drained(ind);
                ch4.weighted = ch4.resample.natural(ind).*0.93 +...
                    ch4.resample.drained(ind).*0.07;
                if ii==1
                    FFI = 100/0.35;
                else
                    FFI = round(100/AAB(j));
                end
                NCB_CC(k,j,i,ii) = C_balance_model(NEE.weighted,NEE_burned(ind),C_loss_fire,...
                    FFI,years_burn_lag(ind), years_2_recovery(ind)) + ch4.weighted;
            end
        end
    end
end

% Increase burn rate
tmp = squeeze(mean(NCB_CC(:,1:30,1,2),2)).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_CC_2050 = [mean(tmp) std(tmp)];

% Increase consumption
tmp = squeeze(mean(NCB_CC(:,1:30,2,1),2)).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_CC_2050(2,:) = [mean(tmp) std(tmp)];

% Increase of both
tmp = squeeze(mean(NCB_CC(:,1:30,2,2),2)).*peat_area./(1e6)./(1e09) .* 30;
NCB_GtC_CC_2050(3,:) = [mean(tmp) std(tmp)];


tmp = squeeze(mean(NCB_CC(:,1:80,1,2),2)).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_CC_2100 = [mean(tmp) std(tmp)];

tmp = squeeze(mean(NCB_CC(:,1:80,2,1),2)).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_CC_2100(2,:) = [mean(tmp) std(tmp)];

tmp = squeeze(mean(NCB_CC(:,1:80,2,2),2)).*peat_area./(1e6)./(1e09) .* 80;
NCB_GtC_CC_2100(3,:) = [mean(tmp) std(tmp)];


x = 1:6;
y = [NCB_GtC_2050(:,1); NCB_GtC_CC_2050(:,1)];

colorPallette = {'#7f7f7f','#ddcc77','#0072b2','#ffde3e','#fc921f','#ed5151'};
figure('color','white','units','inches','position',[2 2 6.5 3])
hold on; box on
for i=1:6
    bar(x(i),y(i),'edgecolor','k','facecolor',colorPallette{i})
end

%--> add error bars (SD)
eh = errorbar(x,y,[NCB_GtC_2050(:,2); NCB_GtC_CC_2050(:,2)]);    
eh.Color = [0 0 0];                            
eh.LineStyle = 'none';


x = 8:13;
y = [NCB_GtC_2100(:,1); NCB_GtC_CC_2100(:,1)];
for i=1:6
    bar(x(i),y(i),'edgecolor','k','facecolor',colorPallette{i})
end

%--> add error bars (SD)
eh = errorbar(x,y,[NCB_GtC_2100(:,2); NCB_GtC_CC_2100(:,2)]);
eh.Color = [0 0 0];
eh.LineStyle = 'none';
eh.Color = [0 0 0];

ylabel('Cumulative NEE + CH_4 (GtC)')

set(gca,'yminortick','on','xtick',[1:7 9:15],'xticklabel','','ylim',[-9 2],'xlim',[0 14])
text(2.5,2.5,'2050','fontweight','bold')
text(10.5,2.5,'2100','fontweight','bold')
text(12.4,1.5,'source','fontweight','bold')
text(12.9,-8.25,'sink','fontweight','bold')
%--> Seperator
plot([7 7],[-12 2.5],'k-','linewidth',1)

xlabels = {'No burn','Current state','Restoration','Increased burn rate',...
    'Increased severity','Full climate change'};
legend(xlabels,'location','southwest')

set(gca,'ylim',[-9 2])
legend(xlabels,'location','southwest')


%% Analysis of FiredPy data and Hugelius histosol fraction
%--> Uses hexscatter.m
%--> https://www.mathworks.com/matlabcentral/fileexchange/45639-hexscatter-m

% Using current directory, so assumes the CSV is in the same folder as the 
%   'Wilkinson_et_al_NCC.m' file.
folder = cd;
file = 'Boreal_firePct_histosolFrac_ecozone - NP.csv';

data = readtable(fullfile(folder,file));

% Combine data on  percent burn from different boreal/temperate
%   non-permafrost regions.
fire_pct = sum([data.fire_pct, data.fired_eurasia_4326_pc,...
    data.fired_russia_china_japan_4326_pc],2,'omitnan');

[histosol_bin,histosol_binEdges] = discretize(data.x_mean,12);
n = histcounts(data.x_mean,histosol_binEdges);
figure('color','white')
bh = boxplot(fire_pct,histosol_bin,'symbol','.','jitter',0.25);
fixBoxplot(bh)
set(gca,'xticklabel',mean([histosol_binEdges(1:end-1)' histosol_binEdges(2:end)'],2),...
    'yscale','log','ylim',[1e-03 1e02])
ylabel('% area burned per 20.75 years')
xlabel('Binned histosol area (%)')

for i=1:12
    if n(i)>99
        offset = 0.2;
    elseif n(i)>999
        offset = 0.3;
    else
        offset = 0.1;
    end
    text(i-2.*offset,150,num2str(n(i)))
end

idx.fire = fire_pct>0;
[fire_bin,fire_binEdges] = discretize(fire_pct(idx.fire),12);
n = histcounts(fire_pct(idx.fire),fire_binEdges);
figure('color','white')
bh = boxplot(data.x_mean(idx.fire),fire_bin,'symbol','.','jitter',0.25);
fixBoxplot(bh)
set(gca,'xticklabel',mean([fire_binEdges(1:end-1)' fire_binEdges(2:end)'],2),...
    'yscale','log','ylim',[0.01 100])
ylabel('Histosol area in regions with fire (%)')
xlabel('Binned fire area per 20.75 years (%)')
for i=1:12
    if n(i)>99
        offset = 0.2;
    elseif n(i)>999
        offset = 0.3;
    else
        offset = 0.1;
    end
    text(i-2.*offset,150,num2str(n(i)))
end

% Mean for ecoregion taking into account area of partial hexagons
ecoregion = unique(data.ECO_NAME);
fire_pct = nan(length(ecoregion),2);
histosol_pct = nan(length(ecoregion),1);
for i=1:length(ecoregion)
    idx.ecoregion = strcmp(data.ECO_NAME,ecoregion{i});
    area_total = sum(data.area(idx.ecoregion),'omitnan');
    W = data.area(idx.ecoregion)./area_total;
    
    fire_pct(i,1) = sum(W.*data.fire_pct(idx.ecoregion),'omitnan');
    fire_pct(i,2) = area_total;
    histosol_pct(i,1) = sum(W.*data.x_mean(idx.ecoregion),'omitnan');
end

T = table(ecoregion,fire_pct,histosol_pct);
writetable(T,fullfile(folder,'Ecoregion_fire_histosol_pct.csv'))


figure('color','white')
hold on; box on
scatter(histosol_pct, fire_pct./20.75,'ko')
hold on; box on
ylabel('Ecozone annual area burned (%)')
xlabel('Ecozone mean histosol cover (%)')
set(gca,'ylim',[-0.05 2],'xminortick','on','yminortick','on')
str = {'Southern Hudson', 'Bay taiga'};
annotation('textbox',[0.7 0.05 0.2 0.3],'String',str,'FitBoxToText','on',...
    'linestyle','none','horizontalalignment','center');

figure('color','white')
hold on; box on
scatter(histosol_pct, 100./(fire_pct./20.75),'ko')
set(gca,'yscale','log','ylim',[50 50000])
ylabel('Fire return interval')
xlabel('Ecozone mean histosol cover (%)')
annotation('textbox',[0.7 0.1 0.2 0.3],'String',str,'FitBoxToText','on',...
    'linestyle','none','horizontalalignment','center');


% Overall fire return interval
%--> Combine data on  percent burn from different boreal/temperate
%       non-permafrost regions.
fire_pct = sum([data.fire_pct, data.fired_eurasia_4326_pc,...
    data.fired_russia_china_japan_4326_pc],2,'omitnan');
W = data.area./sum(data.area,'omitnan'); % weighting
mean_fire_pct = sum(fire_pct.*W,'omitnan');
FRI = 100./(mean_fire_pct./19.75);
mean_fire_pct = sum(fire_pct(data.x_mean>=1).*W(data.x_mean>=1),'omitnan');
FRI_histosol = 100./(mean_fire_pct./19.75);


% Fire return interval by continent/biome
files{1} = 'Boreal_firePct_histosolFrac - Boreal_NA - NP.csv';
files{2} = 'Boreal_firePct_histosolFrac - Boreal_Asia - NP.csv';
files{3} = 'Boreal_firePct_histosolFrac - Boreal_Europe - NP.csv';
files{4} = 'Boreal_firePct_histosolFrac - Temperate_NA - NP.csv';
files{5} = 'Boreal_firePct_histosolFrac - Temperate_Asia - NP.csv';
files{6} = 'Boreal_firePct_histosolFrac - Temperate_Europe - NP.csv';

FRI_biome_continent = nan(6,4);
for i=1:6
    data = readtable(fullfile(folder,files{i}));
    
    % Area weighting
    W = data.area./sum(data.area,'omitnan'); % weighting
    
    % Aggregate fire percent
    fire_pct = sum([data.fire_pct, data.fired_eurasia_4326_pc,...
        data.fired_russia_china_japan_4326_pc],2,'omitnan');
    
    mean_fire_pct = sum(fire_pct.*W,'omitnan');
    FRI_biome_continent(i,1) = 100./(mean_fire_pct./19.75);
    
    mean_fire_pct = sum(fire_pct(data.x_mean>=1).*W(data.x_mean>=1),'omitnan');
    FRI_biome_continent(i,2) = 100./(mean_fire_pct./19.75);
    
    FRI_biome_continent(i,3) = sum(fire_pct.*data.area,'omitnan')./100;
    FRI_biome_continent(i,4) = sum(data.area,'omitnan');
end


file = 'Boreal_firePct_histosolFrac_ecozone_CMI - NP.csv';
data = readtable(fullfile(folder,file));
fire_pct = sum([data.fire_pct, data.fired_eurasia_4326_pc,...
    data.fired_russia_china_japan_4326_pc],2,'omitnan');


figure('color','white')
y = log10(fire_pct);
idx = ~isinf(y) & ~isnan(y) & ~isnan(data.x_CMI_mean);
y = y(idx);
x = data.x_CMI_mean(idx);
hexscatter(x, y, 'res',25)
hold on; box on
ch = colorbar;
set(ch,'TickLabels',int16(10.^get(ch,'Ticks')))
set(gca,'fontsize',12)
xlabel('CMI')
ylabel('log10 fire area (%)')

p = polyfit(x,y,1);
plot([-1 1],p(1).*[-1 1]+p(2),'k--')
text(-0.9,-3.5,char(['r^2 = ' num2str(corr(x,y)^2)]),'fontsize',12)


%--> Relationship between histosol and fire
figure('color','white')
y = log10(fire_pct);
x = log10(data.x_mean);
idx = ~isinf(y) & ~isnan(y) & ~isnan(x) & ~isinf(x);
y = y(idx);
x = x(idx);
hexscatter(x, y, 'res',25)
hold on; box on
ch = colorbar;
set(ch,'TickLabels',int16(10.^get(ch,'Ticks')))
set(gca,'fontsize',12)
xlabel('log10 histosol area (%)')
ylabel('log10 fire area (%)')
text(-2.8,-3.5,char(['r^2 = ' num2str(corr(x,y)^2)]),'fontsize',12)

%--> Relationship between CMI and histosol
figure('color','white')
y = log10(data.x_mean);
x = data.x_CMI_mean;
idx = ~isinf(y) & ~isnan(y) & ~isnan(x) & ~isinf(x);
y = y(idx);
x = x(idx);
hexscatter(x, y, 'res',25)
hold on; box on
ch = colorbar;
set(ch,'TickLabels',int16(10.^get(ch,'Ticks')))
set(gca,'fontsize',12)
ylabel('log10 histosol area (%)')
xlabel('CMI')
p = polyfit(x,y,1);
plot([-1 1],p(1).*[-1 1]+p(2),'k--')
set(gca,'xlim',[-1 1])
text(-0.9,-2.25,char(['r^2 = ' num2str(corr(x,y)^2)]),'fontsize',12)

