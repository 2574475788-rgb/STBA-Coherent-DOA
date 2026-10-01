function cfg = config_v0()
%CONFIG_V0 Fixed configuration for the first ML-CFR validation version.
%
% No toolboxes are required beyond base MATLAB.
%
% Main array:
%   10-sensor two-level nested array
%   p = [0 1 2 3 4 5 11 17 23 29], in units of d=lambda/2.
%
% This first version deliberately studies only K=2 fully coherent sources.
% The purpose is to validate the multi-lag common-factor front end before
% adding weighting, array interpolation, unknown-K processing, etc.

cfg.p = [0 1 2 3 4 5 11 17 23 29];
cfg.M = numel(cfg.p);
cfg.K = 2;

cfg.angleMin = -60;        % degrees
cfg.angleMax =  60;        % degrees
cfg.randomMinSep = 5;      % degrees, for main random-angle experiments

cfg.lags = 1:4;            % pre-registered Version-0 lag set
cfg.omega0 = 0.2*pi;       % normalized digital radian frequency

% Diagnostic DOA back end:
% 1-degree global pair search + continuous variable-projection refinement.
cfg.coarseGrid = cfg.angleMin:1:cfg.angleMax;
cfg.coarseMinSep = 1;      % avoid identical/near-identical grid columns
cfg.refineMinSep = 0.05;   % numerical singularity guard only
cfg.localMaxIter = 400;
cfg.localTolX = 1e-8;
cfg.localTolFun = 1e-10;

cfg.seed = 20260923;

% Numerical tolerances
cfg.epsPhase = 1e-12;
cfg.epsEig = 1e-12;
end
