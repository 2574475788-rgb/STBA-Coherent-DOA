function [thetaHat,info] = fbss_music_ula_estimate(X,p,opts)
%FBSS_MUSIC_ULA_ESTIMATE Forward/backward spatial-smoothing MUSIC baseline.
%
% The source paper describes FBSS but does not state a numerical subarray
% length for its simulations.  This project uses the pre-specified
% maximum-aperture choice L=M-K+1, which leaves Q=K overlapping forward
% subarrays and satisfies Q>=K for decorrelation.  For M=10,K=2, L=9.
% This is an explicit favorable adaptation, not an author-code reproduction.

if nargin<3 || isempty(opts), opts=struct(); end
[M,N]=size(X); p=p(:).';
assert(numel(p)==M,'p must contain one position per row of X.');
assert(max(abs(diff(p)-1))<1e-12, ...
    'FBSS baseline requires a unit-spaced half-wavelength ULA.');
if ~isfield(opts,'K'), opts.K=2; end
if ~isfield(opts,'thetaGrid'), opts.thetaGrid=-60:0.1:60; end
K=opts.K;
if ~isfield(opts,'subarrayLength') || isempty(opts.subarrayLength)
    opts.subarrayLength=M-K+1;
end
L=opts.subarrayLength; Q=M-L+1;
assert(L>K && L<=M,'Require K < subarrayLength <= M.');
assert(Q>=K,'Require at least K translated subarrays for K coherent sources.');

Rf=zeros(L,L);
for q=1:Q
    Xq=X(q:q+L-1,:);
    Rf=Rf+(Xq*Xq')/N;
end
Rf=Rf/Q;
J=fliplr(eye(L));
Rb=J*conj(Rf)*J;
Rfb=(Rf+Rb)/2; Rfb=(Rfb+Rfb')/2;

[thetaHat,Pmusic,musicInfo]=tmr_music_from_cov(Rfb,K,opts.thetaGrid);
thetaHat=sort(thetaHat(:).');

info.subarrayLength=L; info.numSubarrays=Q;
info.Rforward=Rf; info.Rbackward=Rb; info.Rfb=Rfb;
info.musicSpectrum=Pmusic; info.music=musicInfo;
info.reproductionStatus=['Equation-based FBSS adaptation; source paper does not ' ...
    'specify its simulation subarray length.'];
info.subarrayRationale='L=M-K+1 maximizes aperture subject to Q=M-L+1 >= K.';
end
