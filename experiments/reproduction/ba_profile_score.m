function [score, info] = ba_profile_score(R0, Rlags, Bcand)
%BA_PROFILE_SCORE Profiled residual-evidence score for a candidate B.
%
% The SAME working model used by BA-FGLS is used here:
%
%   P z0 = c0 P q + e0
%   zj    = cj q     + ej,   j=1,...,J
%
% with scalar-isotropic residual covariance within each lag channel.
%
% For a fixed candidate direction q, the nuisance complex coefficient c_j
% and scalar residual variance nu_j have closed-form ML estimates. After
% profiling them out, the negative log-likelihood differs only by
%
%   S(q) = d0*log(nu0) + sum_j d*log(nu_j),
%
% where
%   d0 = M^2-2  (projected zero-lag residual DOF after fitting c0)
%   d  = M^2-1  (positive-lag residual DOF after fitting cj).
%
% Lower score is better.
%
% This score has NO empirical threshold and NO tuning parameter. It is a
% model-selection criterion under the same scalar-residual approximation
% already used in BA-FGLS.

[M,M2,J] = size(Rlags);
assert(M==M2,'Rlags must contain square matrices.');
assert(all(size(R0)==[M M]),'R0 size mismatch.');

Pdim = M*M;

q = reshape(Bcand,[],1);
q = q/max(norm(q),eps);

z0 = reshape(R0,[],1);
Z = complex(zeros(Pdim,J));
for jj=1:J
    Z(:,jj)=reshape(Rlags(:,:,jj),[],1);
end

ivec=reshape(eye(M),[],1);
inorm2=real(ivec'*ivec);
projP=@(v) v-ivec*((ivec'*v)/inorm2);

% Positive-lag profiled residuals.
c = q'*Z;
E = Z-q*c;
rss = sum(abs(E).^2,1);
d = max(Pdim-1,1);
nu = rss/d;

% Projected zero-lag profiled residual.
y0=projP(z0);
pq=projP(q);
den0=real(pq'*pq);
if den0<=eps
    c0=0;
else
    c0=(pq'*y0)/den0;
end
e0=y0-c0*pq;
d0=max(Pdim-2,1);
rss0=real(e0'*e0);
nu0=rss0/d0;

% Numerical floor only.
allNu=[nu0,nu];
pos=allNu(allNu>0 & isfinite(allNu));
if isempty(pos)
    floorv=eps;
else
    floorv=max(eps,1e-12*median(pos));
end
nu0=max(nu0,floorv);
nu=max(nu,floorv);

score=d0*log(nu0)+d*sum(log(nu));

info.c0=c0;
info.c=c;
info.nu0=nu0;
info.nu=nu;
info.rss0=rss0;
info.rss=rss;
info.score=score;
end
