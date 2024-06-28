function C_balance = C_balance_model(NEE_prefire,NEE_postfire,...
    C_loss_fire,FFI,years_burn_lag, years_2_recovery)
%==========================================================================
% Simple book-keeping model for estimating peatland NEE over a fire return
%   interval.
%
% NEE_prefire	    - Annual average NEE for a peatland pre-fire; Units: g C m^-2 yr^-1
% NEE_postfire      - Annual average NEE for a peatland post-fire; Units: g C m^-2 yr^-1
% C_loss_fire       - Peatland peat carbon loss due to wildfire; Units: kg C
% FFI               - Fire free interval
% years_burn_lag    - Delay in start of NEE recovery post-fire
% year_2_recovery   - Length of time post-fire for NEE to recover to a pre-fire state.
%
%==========================================================================

if nargin==1
    model_inputs = NEE_prefire;
    NEE_postfire = model_inputs(2);
    C_loss_fire = model_inputs(3);
    FFI = model_inputs(4);
    years_burn_lag = model_inputs(5);
    years_2_recovery = model_inputs(6);
    NEE_prefire = model_inputs(1);
end

% Force values to be integers
FFI = round(FFI);
years_burn_lag = round(years_burn_lag);
years_burn_lag(years_burn_lag<=0) = 1;
years_2_recovery = round(years_2_recovery);

% Initialize variable for annual NEE over a fire free interval (FFI)
NEE_annual = nan(FFI,1);
%--> NEE takes on postfire value during the lag in recovery period
NEE_annual(1:years_burn_lag,1) = NEE_postfire;
%--> Once fully recovered, NEE takes on prefire value
NEE_annual(years_2_recovery:end,1) = NEE_prefire;
%--> Length of time over which recovery takes place
delta_x = (years_2_recovery - years_burn_lag) + 1;
%--> Difference between postfire and prefire NEE values
delta_y = NEE_prefire - NEE_postfire;
%--> NEE recovers linearly between postfire and prefire value
NEE_annual(years_burn_lag:years_2_recovery,1) =  delta_y.*(1:delta_x)./delta_x +...
    NEE_postfire;

% Carbon balance in g C m^-2 yr^-1: positive is source
%--> Average NEE_annual over the fire free interval adjusted for the peat
%       carbon loss from fire.
C_balance = (sum(NEE_annual) + 1000.*C_loss_fire)./FFI;