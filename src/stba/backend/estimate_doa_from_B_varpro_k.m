function [thetaHat, info] = estimate_doa_from_B_varpro_k(B, cfg, K, L, beamWidth, minStartDist)
%ESTIMATE_DOA_FROM_B_VARPRO_K General K-source continuous decoder.
%
% 1) Extract principal coherent spatial vector b from rank-one PSD B.
% 2) Generate multiple deterministic coarse starts with beam search.
% 3) Refine every start by K-dimensional separable nonlinear LS.
% 4) Select the valid refined solution with the smallest profiled residual.
%
% K is treated as known in Version 2A on purpose. Model-order estimation is
% a separate problem and is not mixed into backend qualification.

if nargin < 4 || isempty(L), L=8; end
if nargin < 5 || isempty(beamWidth), beamWidth=80; end
if nargin < 6 || isempty(minStartDist), minStartDist=1.5; end

b=extract_b_from_B(B);

[theta0,cinfo]=coarse_k_beam_candidates( ...
    b,cfg,K,L,beamWidth,minStartDist);

opts=optimset('Display','off', ...
    'MaxIter',cfg.localMaxIter, ...
    'TolX',cfg.localTolX, ...
    'TolFun',cfg.localTolFun);

obj=@(th) varpro_objective_k(th,b,cfg,K);

nCand=size(theta0,1);
thetaRefined=nan(nCand,K);
fvalRefined=inf(nCand,1);
exitflag=nan(nCand,1);
usedFallback=false(nCand,1);

for cc=1:nCand
    th0=theta0(cc,:);
    f0=obj(th0);

    [th,fv,ef]=fminsearch(obj,th0,opts);
    th=sort(th(:).');

    valid = numel(th)==K && ...
        all(th>=cfg.angleMin) && all(th<=cfg.angleMax) && ...
        (K==1 || min(diff(th))>=cfg.refineMinSep) && ...
        isfinite(fv);

    if valid && fv <= f0
        thetaRefined(cc,:)=th;
        fvalRefined(cc)=fv;
    else
        thetaRefined(cc,:)=th0;
        fvalRefined(cc)=f0;
        usedFallback(cc)=true;
    end

    exitflag(cc)=ef;
end

[bestF,bestIdx]=min(fvalRefined);
thetaHat=thetaRefined(bestIdx,:);

% Estimate coherent amplitudes at the selected DOAs for diagnostics.
Ahat=steering_matrix(cfg.p,thetaHat);
betaHat=Ahat\b;
r=b-Ahat*betaHat;

info.b=b;
info.thetaCoarse=theta0;
info.thetaRefined=thetaRefined;
info.fvalRefined=fvalRefined;
info.exitflag=exitflag;
info.usedFallback=usedFallback;
info.bestIndex=bestIdx;
info.bestFval=bestF;
info.betaHat=betaHat;
info.finalResidual=real(r'*r)/max(real(b'*b),eps);
info.coarseInfo=cinfo;
info.K=K;
info.L=L;
info.beamWidth=beamWidth;
info.minStartDist=minStartDist;
end
