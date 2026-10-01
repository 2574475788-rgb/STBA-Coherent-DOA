function [theta2,info] = tmr_vai_fully_augmented_estimate( ...
    X,p,K,sourceGrid,musicGrid,duplicateLagMode)
%TMR_VAI_FULLY_AUGMENTED_ESTIMATE
% Corrected fully-augmented specialization of Zhang et al. TMR-VAI.
%
% Source-paper rule:
% For a hole-free difference coarray, after decorrelation one may directly
% use coarray MUSIC instead of virtual-hole interpolation / Eq. (27).
%
% Here "coarray MUSIC" is implemented with the spatial-smoothing
% construction associated with the paper's cited nested-array reference:
%
%   filled lag vector -> overlapping virtual ULA subarrays
%                     -> spatially smoothed covariance
%                     -> known-K MUSIC.
%
% The two-stage TMR source-covariance refinement is retained:
%   Eq.(17) -> remove off-diagonal source terms -> coarray SS-MUSIC
%   -> theta^(1) -> Eq.(28)/(33) -> repeat -> theta^(2).
%
% No CVX is required because there are no virtual holes to interpolate.

if nargin < 6 || isempty(duplicateLagMode)
    duplicateLagMode='mean';
end

L=size(X,2);
Rx=(X*X')/L;
Rx=(Rx+Rx')/2;

%% Pass 1: Eq. (17)
s1=tmr_source_cov_initial(Rx,p,sourceGrid);

Rdec1=tmr_remove_source_cross_terms(Rx,s1.Agrid,s1.Rs);
co1=tmr_coarray_observation(Rdec1,p,duplicateLagMode);

if ~isempty(co1.holes)
    error(['tmr_vai_fully_augmented_estimate requires a hole-free ' ...
           'difference coarray.']);
end

[Rss1,ss1]=tmr_coarray_ss_covariance(co1);
[theta1,P1,mus1]=tmr_music_from_cov(Rss1,K,musicGrid);

%% Pass 2: Eqs. (28)/(33), overdetermined case
Ahat=tmr_steering(p,theta1);

if K <= numel(p)
    Adag2=(Ahat'*Ahat)\Ahat';
else
    Adag2=Ahat'/(Ahat*Ahat');
end

Rs2=Adag2*Rx*Adag2';
Rs2=(Rs2+Rs2')/2;

Rdec2=tmr_remove_source_cross_terms(Rx,Ahat,Rs2);
co2=tmr_coarray_observation(Rdec2,p,duplicateLagMode);

if ~isempty(co2.holes)
    error('Unexpected holes appeared in the fixed array geometry.');
end

[Rss2,ss2]=tmr_coarray_ss_covariance(co2);
[theta2,P2,mus2]=tmr_music_from_cov(Rss2,K,musicGrid);

info.Rx=Rx;
info.initial=s1;

info.Rdec1=Rdec1;
info.coarray1=co1;
info.Rss1=Rss1;
info.ss1=ss1;
info.theta1=theta1;
info.P1=P1;
info.music1=mus1;

info.Ahat=Ahat;
info.Rs2=Rs2;
info.Rdec2=Rdec2;
info.coarray2=co2;
info.Rss2=Rss2;
info.ss2=ss2;
info.theta2=theta2;
info.P2=P2;
info.music2=mus2;

info.minEig1=ss1.minEig;
info.minEig2=ss2.minEig;
end
