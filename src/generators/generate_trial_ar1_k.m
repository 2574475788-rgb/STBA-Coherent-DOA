function [X, truth] = generate_trial_ar1_k(cfg, N, snrDb, rho, thetaDeg, seed)
%GENERATE_TRIAL_AR1_K K-source fully coherent AR(1)-envelope generator.
%
% Generalizes generate_trial_ar1.m from K=2 to arbitrary cfg.K.
%
% Fully coherent model:
%   s_k[n] = beta_k s[n]
%   X = A(theta) beta s + W
%
% with equal-magnitude coherent coefficients
%   beta_1 = 1,
%   beta_k = exp(j phi_k), k=2,...,K.
%
% Thus the equivalent spatial vector remains
%   b = A(theta) beta
% and all nonzero-lag population matrices are scalar multiples of
%   B = b b^H,
% irrespective of K.
%
% The finite-record array-input SNR is fixed exactly, as in V1A-V1H.

if nargin < 5
    thetaDeg = [];
end
if nargin < 6
    seed = cfg.seed;
end

K = cfg.K;

assert(K >= 1 && K < cfg.M, 'Require 1 <= K < M.');
assert(rho >= 0 && rho <= 1, 'rho must lie in [0,1].');
assert(N > max(cfg.lags), 'N must exceed the maximum lag.');

rng(seed,'twister');

% ----- angles -----
if isempty(thetaDeg)
    ok = false;
    while ~ok
        th = sort(cfg.angleMin + (cfg.angleMax-cfg.angleMin)*rand(1,K));
        if K == 1
            ok = true;
        else
            ok = min(diff(th)) >= cfg.randomMinSep;
        end
    end
    thetaDeg = th;
else
    thetaDeg = sort(thetaDeg(:).');
    assert(numel(thetaDeg)==K,'thetaDeg must contain cfg.K angles.');
end

% ----- fully coherent spatial coefficients -----
% Draw K-1 phases. For K=2 this preserves the same random-variable count
% as the original generator.
phiRel = 2*pi*rand(max(K-1,0),1);
beta = ones(K,1);
if K > 1
    beta(2:end) = exp(1j*phiRel);
end

A = steering_matrix(cfg.p,thetaDeg);
b = A*beta;

% ----- stochastic narrowband common waveform -----
varphi0 = 2*pi*rand;

u0 = (randn + 1j*randn)/sqrt(2);
w = (randn(1,N) + 1j*randn(1,N))/sqrt(2);

u = complex(zeros(1,N));
u(1) = u0;
alpha = sqrt(max(0,1-rho^2));
for n = 2:N
    u(n) = rho*u(n-1) + alpha*w(n);
end

idx = 0:(N-1);
carrier = exp(1j*(cfg.omega0*idx + varphi0));
s = u .* carrier;

Xsig = b*s;

% ----- exact finite-record array-input SNR -----
Psig = norm(Xsig,'fro')^2/(cfg.M*N);
if isinf(snrDb) && snrDb > 0
    sigma2 = 0;
    W = zeros(cfg.M,N);
else
    sigma2 = Psig*10^(-snrDb/10);
    W0 = (randn(cfg.M,N)+1j*randn(cfg.M,N))/sqrt(2);
    W = sqrt(sigma2)*W0;
end

X = Xsig + W;

Btrue = b*b';
Btrue = (Btrue+Btrue')/2;
Btrue = Btrue/real(trace(Btrue));

% empirical source temporal correlation
lags = cfg.lags;
rEmp = complex(zeros(size(lags)));
pEmp = mean(abs(s).^2);
for jj=1:numel(lags)
    tau=lags(jj);
    rEmp(jj)=mean(s(tau+1:end).*conj(s(1:end-tau)))/pEmp;
end

truth.thetaDeg = thetaDeg;
truth.beta = beta;
truth.phiRel = phiRel;
truth.A = A;
truth.b = b;
truth.Btrue = Btrue;
truth.rho = rho;
truth.s = s;
truth.u = u;
truth.Psig = Psig;
truth.sigma2 = sigma2;
truth.rEmp = rEmp;
truth.rTheory = (rho.^lags).*exp(1j*cfg.omega0*lags);
end
