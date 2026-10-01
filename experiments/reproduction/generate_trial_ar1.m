function [X, truth] = generate_trial_ar1(cfg, N, snrDb, rho, thetaDeg, seed)
%GENERATE_TRIAL_AR1 Fully coherent sources with a stochastic narrowband
%Gauss-Markov (AR(1)) complex envelope.
%
%   x[n] = A(theta)*beta*s[n] + noise[n]
%   beta = [1; exp(j*phi)]
%   s[n] = u[n] * exp(j*(omega0*n + varphi0))
%   u[n] = rho*u[n-1] + sqrt(1-rho^2)*w[n]
%
% For stationary unit-variance u[n]:
%   E{s[n]s*[n-tau]} = rho^tau * exp(j*omega0*tau)
% for tau >= 0 (rho real, 0<=rho<=1).
%
% The actual finite-record array input SNR is fixed trial-by-trial:
%   Psig = ||Xsig||_F^2 / (M*N)
%   sigma_n^2 = Psig * 10^(-SNR/10)
%
% IMPORTANT FOR PAIRED rho SWEEPS:
% This function draws exactly the same number/order of random variables for
% every rho. Calling it with the SAME seed but different rho therefore keeps:
%   theta, coherent phase, carrier phase, AR innovations, and standard-normal
%   noise direction paired across rho values.

if nargin < 5
    thetaDeg = [];
end
if nargin < 6
    seed = cfg.seed;
end

assert(rho >= 0 && rho <= 1, 'rho must lie in [0,1].');
assert(N > max(cfg.lags), 'N must exceed the maximum lag.');

rng(seed,'twister');

% ----- angles -----
if isempty(thetaDeg)
    ok = false;
    while ~ok
        th = sort(cfg.angleMin + (cfg.angleMax-cfg.angleMin)*rand(1,cfg.K));
        ok = min(diff(th)) >= cfg.randomMinSep;
    end
    thetaDeg = th;
else
    thetaDeg = sort(thetaDeg(:).');
    assert(numel(thetaDeg)==cfg.K,'thetaDeg must contain K angles.');
end

% ----- coherent spatial coefficients -----
phi = 2*pi*rand;
beta = [1; exp(1j*phi)];
A = steering_matrix(cfg.p,thetaDeg);
b = A*beta;

% ----- stochastic narrowband common waveform -----
varphi0 = 2*pi*rand;

% Draw initial state and ALL innovations even when rho=1 so RNG pairing is
% preserved across rho.
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

% Exact finite-record array-input SNR.
Psig = norm(Xsig,'fro')^2/(cfg.M*N);
if isinf(snrDb) && snrDb>0
    sigma2 = 0;
    W = zeros(cfg.M,N);
else
    sigma2 = Psig*10^(-snrDb/10);

    % Always draw this matrix after the same preceding RNG operations.
    W0 = (randn(cfg.M,N)+1j*randn(cfg.M,N))/sqrt(2);
    W = sqrt(sigma2)*W0;
end

X = Xsig + W;

Btrue = b*b';
Btrue = Btrue/real(trace(Btrue));

% Empirical normalized temporal correlations for diagnostics.
lags = cfg.lags;
rEmp = complex(zeros(size(lags)));
pEmp = mean(abs(s).^2);
for jj = 1:numel(lags)
    tau = lags(jj);
    rEmp(jj) = mean(s(tau+1:end).*conj(s(1:end-tau)))/pEmp;
end

truth.thetaDeg = thetaDeg;
truth.phi = phi;
truth.beta = beta;
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
