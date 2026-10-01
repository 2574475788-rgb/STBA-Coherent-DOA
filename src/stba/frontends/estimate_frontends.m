function out = estimate_frontends(X, cfg)
%ESTIMATE_FRONTENDS Run all Version-0 common-factor front ends on SAME data.

N = size(X,2);
R0 = X*X'/N;
Rlags = lagged_covariances(X, cfg.lags);

[out.B_ZLPCA, out.info_ZLPCA] = zero_lag_pca(R0);
[out.B_SL,    out.info_SL]    = single_lag_cfr(Rlags(:,:,1), cfg);
[out.B_PA,    out.info_PA]    = phase_aligned_mla(Rlags, cfg);
[out.B_ML,    out.info_ML]    = ml_cfr(Rlags, cfg);

out.R0 = R0;
out.Rlags = Rlags;
end
