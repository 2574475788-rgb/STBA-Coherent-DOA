function [B, info] = fgls_cfr(Rlags, cfg)
%FGLS_CFR One-step feasible-GLS common-factor recovery.
%
% Lag-domain rank-one regression:
%   z_j = c_j q + e_j.
%
% Assumption used ONLY for the WLS derivation:
%   E[e_j]=0,
%   Cov(e_j) ~= nu_j I within each stacked lag column.
%
% Pilot:
%   q0 = dominant left singular vector of Z.
%   c_j = q0^H z_j.
%
% Residual variance proxy:
%   nu_j = ||z_j - c_j q0||^2 / (M^2-1).
%
% For fixed c_j,nu_j, the GLS minimizer is
%
%   q_F = [sum_j c_j^* z_j / nu_j] /
%         [sum_j |c_j|^2 / nu_j].
%
% We perform ONE feasible update only. This avoids iterative overfitting and
% keeps the first publishable candidate closed-form and interpretable.

[M,M2,J] = size(Rlags);
assert(M==M2,'Rlags must contain square matrices.');

P = M*M;
Z = complex(zeros(P,J));
for jj=1:J
    Z(:,jj) = reshape(Rlags(:,:,jj),[],1);
end

[U,S,~] = svd(Z,'econ');
q0 = U(:,1);
q0 = q0/norm(q0);

cHat = q0' * Z;  % 1 x J
E = Z - q0*cHat;

nuHat = sum(abs(E).^2,1) / max(P-1,1);

% Numerical floor only; it is relative to the observed residual scale and
% does not act as a statistical tuning parameter.
pos = nuHat(nuHat>0 & isfinite(nuHat));
if isempty(pos)
    nuFloor = eps;
else
    nuFloor = max(eps,1e-8*median(pos));
end
nuUse = max(nuHat,nuFloor);

num = complex(zeros(P,1));
den = 0;
for jj=1:J
    num = num + conj(cHat(jj))*Z(:,jj)/nuUse(jj);
    den = den + abs(cHat(jj))^2/nuUse(jj);
end

if den <= eps
    q = q0;
else
    q = num/den;
    q = q/norm(q);
end

Q = reshape(q,M,M);
[B,pinfo] = project_rank1_psd(Q,cfg);

reliability = abs(cHat).^2 ./ nuUse;
if sum(reliability)>0
    reliabilityNorm = reliability/sum(reliability);
else
    reliabilityNorm = ones(1,J)/J;
end

info.Z = Z;
info.pilotSingularValues = diag(S);
info.q0 = q0;
info.cHat = cHat;
info.nuHat = nuHat;
info.nuUse = nuUse;
info.nuFloor = nuFloor;
info.reliability = reliability;
info.reliabilityNorm = reliabilityNorm;
info.projection = pinfo;
end
