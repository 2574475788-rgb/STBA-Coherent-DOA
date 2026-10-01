function [B, info] = ba_fgls_with_pilot(R0, Rlags, q0, cfg)
%BA_FGLS_WITH_PILOT Bias-annihilated augmented-lag FGLS using a supplied pilot.
%
% This is a diagnostic helper for Version 1D.
% q0 must estimate q = vec(B) and is normalized internally.
%
% The estimator is the same as V1C BA-FGLS; ONLY the pilot changes.

[M,M2,J] = size(Rlags);
assert(M==M2,'Rlags must contain square matrices.');
assert(all(size(R0)==[M M]),'R0 size mismatch.');

Pdim = M*M;

q0 = q0(:);
q0 = q0 / max(norm(q0),eps);

z0 = reshape(R0,[],1);
Z = complex(zeros(Pdim,J));
for jj=1:J
    Z(:,jj) = reshape(Rlags(:,:,jj),[],1);
end

ivec = reshape(eye(M),[],1);
inorm2 = real(ivec'*ivec);

projP = @(v) v - ivec*((ivec'*v)/inorm2);
projI = @(v) ivec*((ivec'*v)/inorm2);

% Positive lags.
cHat = q0' * Z;
E = Z - q0*cHat;
nuHat = sum(abs(E).^2,1) / max(Pdim-1,1);

% Projected zero lag.
y0 = projP(z0);
pq0 = projP(q0);
denC0 = real(pq0'*pq0);
if denC0 <= eps
    c0Hat = 0;
else
    c0Hat = (pq0' * y0) / denC0;
end

e0 = y0 - c0Hat*pq0;
nu0Hat = real(e0'*e0) / max(Pdim-2,1);

% Numerical floors only.
allNu = [nu0Hat,nuHat];
pos = allNu(allNu>0 & isfinite(allNu));
if isempty(pos)
    nuFloor = eps;
else
    nuFloor = max(eps,1e-8*median(pos));
end
nu0 = max(nu0Hat,nuFloor);
nu = max(nuHat,nuFloor);

alpha0 = abs(c0Hat)^2/nu0;
alphaPos = abs(cHat).^2 ./ nu;
beta = sum(alphaPos);

rhs = conj(c0Hat)/nu0 * y0;
for jj=1:J
    rhs = rhs + conj(cHat(jj))/nu(jj)*Z(:,jj);
end

if beta <= eps
    q = q0;
    usedFallback = true;
else
    q = projP(rhs)/(alpha0+beta) + projI(rhs)/beta;
    if norm(q)<=eps || any(~isfinite(q))
        q = q0;
        usedFallback = true;
    else
        q = q/norm(q);
        usedFallback = false;
    end
end

Q = reshape(q,M,M);
[B,pinfo] = project_rank1_psd(Q,cfg);

rel = [alpha0,alphaPos];
if sum(rel)>0
    relNorm = rel/sum(rel);
else
    relNorm = ones(1,J+1)/(J+1);
end

info.c0Hat = c0Hat;
info.cHat = cHat;
info.nu0Hat = nu0Hat;
info.nuHat = nuHat;
info.nu0Use = nu0;
info.nuUse = nu;
info.reliability = rel;
info.reliabilityNorm = relNorm;
info.usedFallback = usedFallback;
info.projection = pinfo;
end
