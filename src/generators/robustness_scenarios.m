function sc = robustness_scenarios()
%ROBUSTNESS_SCENARIOS Pre-registered V3 model-mismatch scenarios.
%
% All scenarios use the same K=2 geometry distribution, beta phases,
% innovation pool and raw noise realization for a given Monte-Carlo seed.
%
% Temporal models:
%   AR1: in-model envelope with rho=0.8.
%   FIR: 7-tap filtered proper-complex Gaussian envelope.  The taps are
%        fixed a priori and normalized to unit energy.  Its lag-1
%        correlation magnitude (~0.815) is close to the AR1 rho=0.8
%        reference, while higher-lag correlations are deliberately
%        non-geometric.
%
% Spatial noise:
%   C_ij = gamma^(|p_i-p_j|), gamma in {0,0.3,0.6}.
%   gamma=0 gives spatially white noise exactly.

h = [1 0.95 0.80 0.60 0.40 0.25 0.12];
h = h / norm(h);

sc(1) = struct('name','AR1-WHITE', ...
    'temporalType','AR1','rho',0.8,'h',[], ...
    'spatialGamma',0.0);

sc(2) = struct('name','AR1-C03', ...
    'temporalType','AR1','rho',0.8,'h',[], ...
    'spatialGamma',0.3);

sc(3) = struct('name','AR1-C06', ...
    'temporalType','AR1','rho',0.8,'h',[], ...
    'spatialGamma',0.6);

sc(4) = struct('name','FIR-WHITE', ...
    'temporalType','FIR','rho',NaN,'h',h, ...
    'spatialGamma',0.0);

sc(5) = struct('name','FIR-C03', ...
    'temporalType','FIR','rho',NaN,'h',h, ...
    'spatialGamma',0.3);

sc(6) = struct('name','FIR-C06', ...
    'temporalType','FIR','rho',NaN,'h',h, ...
    'spatialGamma',0.6);
end
