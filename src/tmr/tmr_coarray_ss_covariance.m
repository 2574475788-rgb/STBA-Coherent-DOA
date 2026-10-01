function [Rss,info] = tmr_coarray_ss_covariance(co)
%TMR_COARRAY_SS_COVARIANCE Spatially-smoothed coarray covariance.
%
% Input co.lagVals must correspond to contiguous lags -P:P.
%
% For a filled difference coarray, coarray MUSIC in the Pal-Vaidyanathan
% construction uses overlapping virtual ULA subarrays and averages their
% outer products:
%
%   Rss = (1/(P+1)) sum_q z_q z_q^H
%
% where each z_q contains P+1 consecutive virtual-lag samples.
%
% This guarantees a Hermitian PSD matrix (up to roundoff) and performs the
% rank-restoration step used by spatially smoothed coarray MUSIC.
%
% IMPORTANT:
% Do not replace this by eigendecomposition of the raw Toeplitz matrix
% formed directly from the lag vector. Under finite samples / imperfect
% decorrelation that raw matrix can be strongly indefinite.

lags = co.Vlags(:).';
z = co.lagVals(:);

P = max(abs(lags));

if ~isequal(lags,-P:P)
    error('Spatial smoothing requires a contiguous symmetric lag set -P:P.');
end

Lsub = P+1;
nSub = P+1;

Rss = complex(zeros(Lsub,Lsub));

for q=1:nSub
    zq = z(q:q+Lsub-1);
    Rss = Rss + zq*zq';
end

Rss = Rss/nSub;
Rss = (Rss+Rss')/2;

ev = sort(real(eig(Rss)),'descend');

info.P=P;
info.Lsub=Lsub;
info.nSub=nSub;
info.eigenvalues=ev;
info.minEig=min(ev);
info.traceReal=real(trace(Rss));
end
