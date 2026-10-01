function [B, info] = zero_lag_pca(R0)
%ZERO_LAG_PCA Rank-one signal-factor estimate from zero-lag covariance.
%
% This is a diagnostic control. For the ideal fully coherent model,
% the dominant eigenvector of R0 = b*b^H + sigma^2 I is b.

H = (R0 + R0')/2;
[V,D] = eig(H);
evals = real(diag(D));
[lambda1,idx] = max(evals);
u1 = V(:,idx);

B = u1*u1';
B = B / real(trace(B));

info.lambda1 = lambda1;
info.eigvals = evals;
end
