function Rlags = lagged_covariances(X, lags)
%LAGGED_COVARIANCES Sample nonzero-lag cross-correlation matrices.
%
% R_tau = 1/(N-tau) * X(:,tau+1:N) * X(:,1:N-tau)^H
%
% IMPORTANT:
% Do NOT Hermitian-symmetrize R_tau here. Finite-sample R_tau is generally
% non-Hermitian; preserving this is part of the intended model.

[M,N] = size(X);
J = numel(lags);
Rlags = complex(zeros(M,M,J));

for jj = 1:J
    tau = lags(jj);
    assert(tau >= 1 && tau < N, 'Each lag must satisfy 1 <= tau < N.');
    Rlags(:,:,jj) = X(:,tau+1:N) * X(:,1:N-tau)' / (N-tau);
end
end
