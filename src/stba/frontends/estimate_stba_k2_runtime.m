function [thetaHat,info] = estimate_stba_k2_runtime(X,cfg)
%ESTIMATE_STBA_K2_RUNTIME Minimal frozen ST-BA end-to-end K=2 estimator.
%
% This wrapper is ONLY for fair runtime measurement.
%
% It executes exactly the frozen ST-BA path:
%   X
%    -> stability_adaptive_ba_fgls
%    -> frozen K=2 VP decoder
%
% It deliberately does NOT call build_frozen_method_bank(), because that
% helper also computes seven unrelated internal baselines and would
% over-count ST-BA runtime.
%
% Frozen K=2 decoder:
%   L=8
%   beamWidth=80
%   minStartDist=1.5 deg
%
% No algorithmic parameter is changed.

[B,stInfo]=stability_adaptive_ba_fgls(X,cfg);

[thetaHat,vpInfo]=estimate_doa_from_B_varpro_k( ...
    B,cfg,2,8,80,1.5);

thetaHat=sort(thetaHat(:).');

info.B=B;
info.stability=stInfo;
info.vp=vpInfo;
end
