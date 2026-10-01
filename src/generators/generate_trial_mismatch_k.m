function [X,truth] = generate_trial_mismatch_k( ...
    cfg,N,snrDb,scenario,thetaDeg,seed)
%GENERATE_TRIAL_MISMATCH_K Frozen-algorithm robustness generator.
%
% Fully coherent model:
%   X = b s + W,  b = A(theta) beta,  B = b b^H.
%
% The algorithm is NOT changed.  Only the data-generating model changes.
%
% Spatial noise:
%   E[w(:,n) w(:,n)^H] = sigma2 * C,
%   C_ij = gamma^(|p_i-p_j|).
% This is a valid principal submatrix of a stationary exponential
% correlation kernel on the integer half-wavelength grid.
%
% Temporal model:
%   AR1 : u[n]=rho u[n-1] + sqrt(1-rho^2)e[n].
%   FIR : u[n]=sum_k h[k] e[n-k], ||h||_2=1.
% Both are modulated by the same narrowband carrier exp(j omega0 n).
%
% Pairing / reproducibility:
% For a fixed seed, geometry, coherent phases, common innovation pool and
% raw spatial-noise draw are shared across all mismatch scenarios.
%
% IMPORTANT exact-regression property:
% For scenario AR1-WHITE with rho=0.8, this function reproduces
% generate_trial_ar1_k exactly for the same seed (up to roundoff).  Extra
% FIR innovations are drawn only AFTER the raw MxN noise matrix, so they do
% not perturb the legacy random draw order.

if nargin < 5, thetaDeg=[]; end
if nargin < 6, seed=cfg.seed; end

K=cfg.K;
M=cfg.M;

assert(N > max(cfg.lags),'N must exceed maximum lag.');
assert(scenario.spatialGamma>=0 && scenario.spatialGamma<1, ...
    'Require 0 <= spatialGamma < 1.');

rng(seed,'twister');

%% Geometry: identical to generate_trial_ar1_k
if isempty(thetaDeg)
    ok=false;
    while ~ok
        th=sort(cfg.angleMin + (cfg.angleMax-cfg.angleMin)*rand(1,K));
        if K==1
            ok=true;
        else
            ok=min(diff(th))>=cfg.randomMinSep;
        end
    end
    thetaDeg=th;
else
    thetaDeg=sort(thetaDeg(:).');
    assert(numel(thetaDeg)==K,'thetaDeg must contain cfg.K angles.');
end

%% Coherent coefficients: identical draw order
phiRel=2*pi*rand(max(K-1,0),1);
beta=ones(K,1);
if K>1
    beta(2:end)=exp(1j*phiRel);
end

A=steering_matrix(cfg.p,thetaDeg);
b=A*beta;

%% Common waveform random pool
varphi0=2*pi*rand;
u0=(randn + 1j*randn)/sqrt(2);

% Legacy AR1 innovations are drawn here.
wCore=(randn(1,N)+1j*randn(1,N))/sqrt(2);

% Legacy raw spatial noise is drawn at the same point in the RNG stream.
W0=(randn(M,N)+1j*randn(M,N))/sqrt(2);

% Extra FIR innovations are drawn only after W0, preserving exact legacy
% AR1-WHITE output for the same seed.
Lh=7;
wExtra=(randn(1,Lh-1)+1j*randn(1,Lh-1))/sqrt(2);

switch upper(scenario.temporalType)
    case 'AR1'
        rho=scenario.rho;
        assert(rho>=0 && rho<=1,'rho must lie in [0,1].');

        u=complex(zeros(1,N));
        u(1)=u0;
        alpha=sqrt(max(0,1-rho^2));
        for n=2:N
            u(n)=rho*u(n-1)+alpha*wCore(n);
        end

        rBase=(rho.^cfg.lags);

    case 'FIR'
        h=scenario.h(:).';
        assert(numel(h)==Lh,'Frozen FIR mismatch uses 7 taps.');
        h=h/norm(h);

        wFull=[wExtra,wCore];
        y=filter(h,1,wFull);
        u=y(Lh:(Lh+N-1));

        rBase=complex(zeros(size(cfg.lags)));
        den=sum(abs(h).^2);
        for jj=1:numel(cfg.lags)
            tau=cfg.lags(jj);
            if tau < numel(h)
                rBase(jj)=sum(h(tau+1:end).*conj(h(1:end-tau)))/den;
            else
                rBase(jj)=0;
            end
        end

    otherwise
        error('Unknown temporalType: %s',scenario.temporalType);
end

idx=0:(N-1);
carrier=exp(1j*(cfg.omega0*idx + varphi0));
s=u.*carrier;
Xsig=b*s;

%% Exact finite-record signal power; target expected array-input SNR
Psig=norm(Xsig,'fro')^2/(M*N);

if isinf(snrDb) && snrDb>0
    sigma2=0;
else
    sigma2=Psig*10^(-snrDb/10);
end

%% Spatial covariance mismatch
gamma=scenario.spatialGamma;
if gamma==0
    C=eye(M);
else
    dist=abs(cfg.p(:)-cfg.p(:).');
    C=gamma.^dist;
    C=(C+C')/2;
end

% Positive definite for gamma<1 on distinct sensor positions.
Lc=chol(C,'lower');

if sigma2==0
    W=zeros(M,N);
else
    W=sqrt(sigma2)*(Lc*W0);
end

X=Xsig+W;

%% Truth
Btrue=b*b';
Btrue=(Btrue+Btrue')/2;
Btrue=Btrue/real(trace(Btrue));

pEmp=mean(abs(s).^2);
rEmp=complex(zeros(size(cfg.lags)));
for jj=1:numel(cfg.lags)
    tau=cfg.lags(jj);
    rEmp(jj)=mean(s(tau+1:end).*conj(s(1:end-tau)))/pEmp;
end

rTheory=rBase.*exp(1j*cfg.omega0*cfg.lags);

truth.thetaDeg=thetaDeg;
truth.beta=beta;
truth.phiRel=phiRel;
truth.A=A;
truth.b=b;
truth.Btrue=Btrue;
truth.s=s;
truth.u=u;
truth.Psig=Psig;
truth.sigma2=sigma2;
truth.Cnoise=C;
truth.spatialGamma=gamma;
truth.temporalType=scenario.temporalType;
truth.rEmp=rEmp;
truth.rTheory=rTheory;

if sigma2>0
    truth.noiseSampleCov=(W*W')/(N*sigma2);
else
    truth.noiseSampleCov=zeros(M);
end

% Small signatures for pairing audits without storing the full random pool.
truth.rawNoiseSignature=[W0(1,1),W0(end,end)];
truth.waveInnovationSignature=[wCore(1),wCore(end),wExtra(1),wExtra(end)];
end
