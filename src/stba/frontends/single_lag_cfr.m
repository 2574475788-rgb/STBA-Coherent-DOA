function [B, info] = single_lag_cfr(R1, cfg)
%SINGLE_LAG_CFR Single-lag structural control.
% Uses exactly the same rank-one PSD projection as ML-CFR.

[B, pinfo] = project_rank1_psd(R1, cfg);
info = pinfo;
end
