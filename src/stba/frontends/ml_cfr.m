function [B, info] = ml_cfr(Rlags, cfg)
%ML_CFR Multi-Lag Common-Factor Recovery.
%
% Stack vec(R_tau) as columns:
%   Z = [vec(R_tau1), ..., vec(R_tauJ)]
%
% Under the ideal fully coherent model:
%   R_tau = r_s(tau) * b*b^H
% hence rank(Z)=1.
%
% The leading left singular vector estimates vec(b*b^H) up to a complex
% scalar. We reshape it and impose the Hermitian PSD rank-one structure.

[M,M2,J] = size(Rlags);
assert(M == M2, 'Rlags must contain square matrices.');

Z = complex(zeros(M*M,J));
for jj = 1:J
    Z(:,jj) = reshape(Rlags(:,:,jj), [], 1);
end

[U,S,V] = svd(Z, 'econ'); %#ok<ASGLU>
svals = diag(S);

q = U(:,1);
Q = reshape(q, M, M);

[B, pinfo] = project_rank1_psd(Q, cfg);

info.Z = Z;
info.singularValues = svals;
if numel(svals) >= 2 && svals(1) > 0
    info.s2OverS1 = svals(2)/svals(1);
else
    info.s2OverS1 = 0;
end
info.projection = pinfo;
end
