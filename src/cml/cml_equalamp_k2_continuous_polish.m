function [thetaHat, info] = cml_equalamp_k2_continuous_polish(X,p,userOpts)
%CML_EQUALAMP_K2_CONTINUOUS_POLISH Continuous local polish of frozen grid/fine CML.
%
% This wrapper does NOT replace cml_equalamp_k2_estimate. It first calls the
% frozen source-faithful grid/fine CML implementation, then locally maximizes
% the SAME concentrated equal-amplitude K=2 likelihood
%   J=(b^H Rhat b)/(b^H b),
%   b=a(theta1)+exp(j*phi)a(theta2),
% from that data-only CML solution.  The physical steering convention used by
% the ML-CFR generator/project is
%   a(theta)=exp(s*j*pi*p*sin(theta)), s in {-1,+1}.
% The sign is inferred data-only from which convention gives the larger
% likelihood at the frozen CML DOA solution; no truth is used.
%
% If the local search is invalid or fails to improve the likelihood, the
% frozen grid/fine solution is retained.  No true DOA is accepted as input.
% Base MATLAB only (fminsearch).

if nargin<3 || isempty(userOpts), userOpts=struct(); end
opts=struct('minSepDeg',0.05,'angleMinDeg',-60,'angleMaxDeg',60, ...
    'maxIter',500,'maxFunEvals',2000,'tolX',1e-9,'tolFun',1e-11, ...
    'phiInitStepDeg',1,'acceptTolRel',1e-11);
fn=fieldnames(userOpts);
for k=1:numel(fn)
    assert(isfield(opts,fn{k}),'Unknown continuous-CML option: %s',fn{k});
    opts.(fn{k})=userOpts.(fn{k});
end
validateattributes(X,{'numeric'},{'2d','nonempty','finite'});
p=double(p(:));
assert(size(X,1)==numel(p),'Rows of X must match numel(p).');

% Frozen source-faithful grid/fine CML. Runtime includes this full stage.
tGrid=tic;
[thetaGrid,gridInfo]=cml_equalamp_k2_estimate(X,p,struct());
gridSec=toc(tGrid);
thetaGrid=sort(double(thetaGrid(:).'));
assert(numel(thetaGrid)==2 && all(isfinite(thetaGrid)), ...
    'Frozen CML returned an invalid K=2 estimate.');

N=size(X,2);
Rhat=(X*X')/N;
Rhat=(Rhat+Rhat')/2;
phiGrid=-180:opts.phiInitStepDeg:(180-opts.phiInitStepDeg);
[signUsed,phi0,scoreGrid,scoreNeg,scorePos]=inferSteeringSign(thetaGrid,Rhat,p,phiGrid);
x0=[thetaGrid phi0];

fmOpts=optimset('Display','off','MaxIter',opts.maxIter, ...
    'MaxFunEvals',opts.maxFunEvals,'TolX',opts.tolX,'TolFun',opts.tolFun);
tPolish=tic;
[xCand,negCand,exitflag,output]=fminsearch(@(x) penalizedObjective(x,Rhat,p,opts,signUsed),x0,fmOpts);
polishSec=toc(tPolish);
thetaCand=sort(double(xCand(1:2)));
phiCand=double(xCand(3));
valid=isfinite(negCand) && all(isfinite(thetaCand)) && isfinite(phiCand) && ...
    thetaCand(1)>=opts.angleMinDeg && thetaCand(2)<=opts.angleMaxDeg && ...
    diff(thetaCand)>=opts.minSepDeg;
if valid
    scoreCand=cmlScore(thetaCand(1),thetaCand(2),phiCand,Rhat,p,signUsed);
else
    scoreCand=-Inf;
end
acceptTol=opts.acceptTolRel*max(1,abs(scoreGrid));
accepted=valid && scoreCand>=scoreGrid-acceptTol;
if accepted
    thetaHat=thetaCand;
    scoreFinal=scoreCand;
    phiFinal=wrap180(phiCand);
else
    thetaHat=thetaGrid;
    scoreFinal=scoreGrid;
    phiFinal=wrap180(phi0);
end

thetaHat=sort(thetaHat(:).');
info=struct();
info.version='Phase-2Q-continuous-CML-polish-v1';
info.baseEstimator='cml_equalamp_k2_estimate';
info.noTruthInput=true;
info.steeringExponentSign=signUsed;
info.steeringConvention='exp(sign*j*pi*p*sin(theta)); sign inferred from frozen CML output, no truth';
info.scoreAtGridNegativeSign=scoreNeg;
info.scoreAtGridPositiveSign=scorePos;
info.thetaGridDeg=thetaGrid;
info.thetaCandidateDeg=thetaCand;
info.thetaFinalDeg=thetaHat;
info.phiInitDeg=wrap180(phi0);
info.phiCandidateDeg=wrap180(phiCand);
info.phiFinalDeg=phiFinal;
info.scoreGrid=scoreGrid;
info.scoreCandidate=scoreCand;
info.scoreFinal=scoreFinal;
info.scoreGain=scoreFinal-scoreGrid;
info.accepted=accepted;
info.exitflag=exitflag;
info.output=output;
info.gridSec=gridSec;
info.polishSec=polishSec;
info.totalSec=gridSec+polishSec;
info.gridInfo=gridInfo;
end

function f=penalizedObjective(x,Rhat,p,opts,signUsed)
if numel(x)~=3 || any(~isfinite(x))
    f=1e12; return;
end
t1=x(1); t2=x(2);
vLow=max(opts.angleMinDeg-t1,0);
vHigh=max(t2-opts.angleMaxDeg,0);
vOrder=max(opts.minSepDeg-(t2-t1),0);
% Keep ordering explicit rather than sorting inside the objective.
if vLow>0 || vHigh>0 || vOrder>0
    f=1e6 + 1e6*(vLow^2+vHigh^2+vOrder^2);
    return;
end
s=cmlScore(t1,t2,x(3),Rhat,p,signUsed);
if ~isfinite(s), f=1e12; else, f=-s; end
end

function [signUsed,phi0,scoreBest,scoreNeg,scorePos]=inferSteeringSign(thetaGrid,Rhat,p,phiGrid)
scoresNeg=zeros(size(phiGrid));scoresPos=zeros(size(phiGrid));
for ip=1:numel(phiGrid)
    scoresNeg(ip)=cmlScore(thetaGrid(1),thetaGrid(2),phiGrid(ip),Rhat,p,-1);
    scoresPos(ip)=cmlScore(thetaGrid(1),thetaGrid(2),phiGrid(ip),Rhat,p,+1);
end
[scoreNeg,iNeg]=max(scoresNeg);[scorePos,iPos]=max(scoresPos);
if scoreNeg>=scorePos
    signUsed=-1;phi0=phiGrid(iNeg);scoreBest=scoreNeg;
else
    signUsed=+1;phi0=phiGrid(iPos);scoreBest=scorePos;
end
end

function s=cmlScore(t1,t2,phiDeg,Rhat,p,signUsed)
a1=exp(signUsed*1j*pi*p*sind(t1));
a2=exp(signUsed*1j*pi*p*sind(t2));
b=a1+exp(1j*pi*phiDeg/180)*a2;
den=real(b'*b);
if ~(isfinite(den) && den>0)
    s=-Inf; return;
end
s=real(b'*Rhat*b)/den;
end

function x=wrap180(x)
x=mod(x+180,360)-180;
end
