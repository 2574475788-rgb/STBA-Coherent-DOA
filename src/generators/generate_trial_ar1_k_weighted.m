function [X, truth] = generate_trial_ar1_k_weighted(cfg,N,snrDb,rho,thetaDeg,seed,betaMag)
%GENERATE_TRIAL_AR1_K_WEIGHTED Fully coherent AR(1) generator with fixed magnitudes.
%
% Model:
%   s_k[n] = beta_k s[n],   beta_k = betaMag(k) exp(j phi_k),
%   beta_1 has zero reference phase and betaMag is supplied by the caller.
%
% Design requirement for the unequal-amplitude experiment:
%   * same seed + same thetaDeg + different betaMag reuses the same relative
%     phases, AR innovations, carrier phase and normalized noise realization;
%   * noise is rescaled per condition to preserve the exact finite-record
%     array-input SNR, so amplitude ratio is not confounded with total SNR;
%   * betaMag=ones(K,1) is bitwise-compatible with generate_trial_ar1_k.

if nargin < 5, thetaDeg=[]; end
if nargin < 6 || isempty(seed), seed=cfg.seed; end
K=cfg.K;
if nargin < 7 || isempty(betaMag), betaMag=ones(K,1); end
betaMag=betaMag(:);
assert(numel(betaMag)==K,'betaMag must contain cfg.K nonnegative magnitudes.');
assert(all(isfinite(betaMag)) && all(betaMag>=0) && any(betaMag>0), ...
    'betaMag must be finite, nonnegative, and not all zero.');
assert(K>=1 && K<cfg.M,'Require 1 <= K < M.');
assert(rho>=0 && rho<=1,'rho must lie in [0,1].');
assert(N>max(cfg.lags),'N must exceed the maximum lag.');

rng(seed,'twister');

% Angles: identical random-number consumption to generate_trial_ar1_k.
if isempty(thetaDeg)
    ok=false;
    while ~ok
        th=sort(cfg.angleMin+(cfg.angleMax-cfg.angleMin)*rand(1,K));
        if K==1, ok=true; else, ok=min(diff(th))>=cfg.randomMinSep; end
    end
    thetaDeg=th;
else
    thetaDeg=sort(thetaDeg(:).');
    assert(numel(thetaDeg)==K,'thetaDeg must contain cfg.K angles.');
end

% Relative phases: identical draw sequence to generate_trial_ar1_k.
phiRel=2*pi*rand(max(K-1,0),1);
beta=complex(betaMag);
if K>1, beta(2:end)=betaMag(2:end).*exp(1j*phiRel); end
beta(1)=betaMag(1);

A=steering_matrix(cfg.p,thetaDeg); b=A*beta;

% Common stochastic narrowband waveform: identical draw sequence.
varphi0=2*pi*rand;
u0=(randn+1j*randn)/sqrt(2);
w=(randn(1,N)+1j*randn(1,N))/sqrt(2);
u=complex(zeros(1,N)); u(1)=u0;
alpha=sqrt(max(0,1-rho^2));
for n=2:N, u(n)=rho*u(n-1)+alpha*w(n); end
idx=0:(N-1); carrier=exp(1j*(cfg.omega0*idx+varphi0)); s=u.*carrier;
Xsig=b*s;

% Exact finite-record array-input SNR; W0 is paired across betaMag values.
Psig=norm(Xsig,'fro')^2/(cfg.M*N);
if isinf(snrDb) && snrDb>0
    sigma2=0; W=zeros(cfg.M,N);
else
    sigma2=Psig*10^(-snrDb/10);
    W0=(randn(cfg.M,N)+1j*randn(cfg.M,N))/sqrt(2);
    W=sqrt(sigma2)*W0;
end
X=Xsig+W;

Btrue=b*b'; Btrue=(Btrue+Btrue')/2; Btrue=Btrue/real(trace(Btrue));
lags=cfg.lags; rEmp=complex(zeros(size(lags))); pEmp=mean(abs(s).^2);
for jj=1:numel(lags)
    tau=lags(jj); rEmp(jj)=mean(s(tau+1:end).*conj(s(1:end-tau)))/pEmp;
end
truth.thetaDeg=thetaDeg; truth.beta=beta; truth.betaMag=betaMag;
truth.phiRel=phiRel; truth.A=A; truth.b=b; truth.Btrue=Btrue;
truth.rho=rho; truth.s=s; truth.u=u; truth.Psig=Psig; truth.sigma2=sigma2;
truth.rEmp=rEmp; truth.rTheory=(rho.^lags).*exp(1j*cfg.omega0*lags);
end
