function out = coherent_equalamp_constrained_crb_k2(cfg,thetaDeg,beta,s,sigma2)
%COHERENT_EQUALAMP_CONSTRAINED_CRB_K2
% Exact-model conditional deterministic CRB for the project's K=2 model.
%
% Model:
%   X(:,n) = b(theta,phi) s(n) + W(:,n)
%   b = a(theta1) + exp(j*phi) a(theta2)
%   beta = [1; exp(j*phi)]
%   W(:,n) ~ CN(0,sigma2 I)
%
% Nuisance parameters:
%   - unknown common deterministic waveform s(n), complex, n=1..N
%   - unknown relative phase phi
%
% Parameters of interest:
%   theta = [theta1,theta2], in radians for differentiation.
%
% After eliminating the unknown complex waveform, the FIM for
% eta=[theta1,theta2,phi] is
%
%   J_eta = 2*||s||^2/sigma2 * Re{ G^H Pi_b_perp G },
%
% where
%   G = [ db/dtheta1, db/dtheta2, db/dphi ],
%   Pi_b_perp = I - b b^H/(b^H b).
%
% The unknown relative phase phi is then removed by a Schur complement.
%
% This bound matches the exact equal-magnitude fully coherent K=2 model
% used by generate_trial_ar1_k. It is distinct from the broader
% deterministic-source CRB that treats the K source waveforms as
% arbitrary nuisance rows.

thetaDeg=thetaDeg(:).';
beta=beta(:);
s=s(:).';

if numel(thetaDeg)~=2 || numel(beta)~=2
    error('coherent_equalamp_constrained_crb_k2:Size', ...
        'Require exactly two DOAs and two coherent coefficients.');
end
if any(diff(thetaDeg)<=0)
    error('coherent_equalamp_constrained_crb_k2:Order', ...
        'thetaDeg must be strictly increasing; formal generator already satisfies this.');
end
if abs(beta(1)-1)>1e-10 || abs(abs(beta(2))-1)>1e-10
    error('coherent_equalamp_constrained_crb_k2:BetaModel', ...
        'Exact-model helper requires beta=[1;exp(j*phi)].');
end
if ~(isfinite(sigma2) && sigma2>0)
    error('coherent_equalamp_constrained_crb_k2:Noise', ...
        'sigma2 must be finite and positive.');
end

p=cfg.p(:);
A=steering_matrix(cfg.p,thetaDeg);
thetaRad=thetaDeg*pi/180;

D=complex(zeros(cfg.M,2));
for k=1:2
    D(:,k)=(-1j*pi*p*cos(thetaRad(k))).*A(:,k);
end

b=A*beta;
bn2=real(b'*b);
Pi=eye(cfg.M)-(b*b')/bn2;
Pi=(Pi+Pi')/2;

G=complex(zeros(cfg.M,3));
G(:,1)=beta(1)*D(:,1);
G(:,2)=beta(2)*D(:,2);
G(:,3)=1j*beta(2)*A(:,2);

Es=sum(abs(s).^2);

Jeta=(2*Es/sigma2)*real(G'*Pi*G);
Jeta=(Jeta+Jeta.')/2;

Jtt=Jeta(1:2,1:2);
Jtp=Jeta(1:2,3);
Jpp=Jeta(3,3);

if ~(isfinite(Jpp) && Jpp>0)
    error('coherent_equalamp_constrained_crb_k2:PhaseFIM', ...
        'Relative-phase information is non-positive or non-finite.');
end

Jtheta=Jtt-(Jtp*Jtp.')/Jpp;
Jtheta=(Jtheta+Jtheta.')/2;

eigJeta=eig(Jeta);
eigJtheta=eig(Jtheta);

minEigJeta=min(real(eigJeta));
minEigJtheta=min(real(eigJtheta));
rcondJeta=rcond(Jeta);
rcondJtheta=rcond(Jtheta);

usedPinv=false;
if ~isfinite(rcondJtheta) || rcondJtheta<1e-14 || minEigJtheta<=0
    C=pinv(Jtheta);
    usedPinv=true;
else
    C=Jtheta\eye(2);
end
C=(C+C.')/2;

out.thetaDeg=thetaDeg;
out.beta=beta;
out.phi=angle(beta(2));
out.sigma2=sigma2;
out.sourceEnergy=Es;

out.A=A;
out.D=D;
out.b=b;
out.PiPerp=Pi;
out.G=G;

out.Jeta=Jeta;
out.Jtheta=Jtheta;
out.Jtt=Jtt;
out.Jtp=Jtp;
out.Jpp=Jpp;

out.Crad2=C;
out.Cdeg2=C*(180/pi)^2;
out.rmseRad=sqrt(max(real(trace(C))/2,0));
out.rmseDeg=out.rmseRad*180/pi;

out.minEigJeta=minEigJeta;
out.minEigJtheta=minEigJtheta;
out.rcondJeta=rcondJeta;
out.rcondJtheta=rcondJtheta;
out.usedPinv=usedPinv;
end
