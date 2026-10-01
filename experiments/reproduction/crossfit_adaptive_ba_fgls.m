function [B, info] = crossfit_adaptive_ba_fgls(X, cfg)
%CROSSFIT_ADAPTIVE_BA_FGLS Two-way cross-fitted pilot selector.
%
% Why:
%   The V1E in-sample profile score reuses the same data both to construct
%   BA-POS / BA-ZLP and to score them. This creates an asymmetric "double
%   dipping" advantage for BA-ZLP because its pilot is extracted from R0 and
%   the score itself contains a zero-lag term.
%
% Remedy:
%   Split snapshots into two non-overlapping contiguous halves A and B.
%
%   Fold A->B:
%       build BA-POS and BA-ZLP using A;
%       score BOTH candidates only on B.
%
%   Fold B->A:
%       build using B;
%       score only on A.
%
%   Sum held-out scores and select the smaller one.
%
%   Finally, recompute the selected candidate on ALL snapshots.
%
% No rho, ground truth, hand threshold, tuned switching point, or extra
% statistical model is introduced.
%
% Note:
%   Signal samples can remain temporally correlated across the split; the
%   key purpose here is to prevent the SAME sample covariance realization
%   from being used simultaneously for fitting and model selection.

N = size(X,2);
nA = floor(N/2);
assert(nA > max(cfg.lags) && (N-nA) > max(cfg.lags), ...
    'Each cross-fit half must contain more samples than max(cfg.lags).');

XA = X(:,1:nA);
XB = X(:,nA+1:end);

% A-trained candidates, B validation.
[candA, outA] = build_candidates(XA,cfg);
outB = estimate_frontends(XB,cfg);

[sPosAB,~] = ba_profile_score(outB.R0,outB.Rlags,candA.Bpos);
[sZlpAB,~] = ba_profile_score(outB.R0,outB.Rlags,candA.Bzlp);

% B-trained candidates, A validation.
[candB, ~] = build_candidates(XB,cfg);

[sPosBA,~] = ba_profile_score(outA.R0,outA.Rlags,candB.Bpos);
[sZlpBA,~] = ba_profile_score(outA.R0,outA.Rlags,candB.Bzlp);

scorePos = sPosAB + sPosBA;
scoreZlp = sZlpAB + sZlpBA;

% Refit both candidates using all snapshots AFTER selection evidence has been
% computed from held-out folds.
[candFull, outFull] = build_candidates(X,cfg);

if scoreZlp < scorePos
    B = candFull.Bzlp;
    selected = 'ZLP';
    selectedIndex = 2;
else
    B = candFull.Bpos;
    selected = 'POS';
    selectedIndex = 1;
end

info.selected = selected;
info.selectedIndex = selectedIndex;
info.scorePos = scorePos;
info.scoreZlp = scoreZlp;
info.scoreMargin = scoreZlp-scorePos; % negative => ZLP
info.foldScores = [sPosAB,sZlpAB; sPosBA,sZlpBA];
info.full.Bpos = candFull.Bpos;
info.full.Bzlp = candFull.Bzlp;
info.full.pos = candFull.ipos;
info.full.zlp = candFull.izlp;
info.outFull = outFull;
end

function [cand,out] = build_candidates(X,cfg)
out = estimate_frontends(X,cfg);

[M,~,J] = size(out.Rlags);
Z = complex(zeros(M*M,J));
for jj=1:J
    Z(:,jj) = reshape(out.Rlags(:,:,jj),[],1);
end

% Positive-lag SVD pilot.
[U,~,~] = svd(Z,'econ');
qPos = U(:,1);
qPos = qPos/max(norm(qPos),eps);

% Zero-lag PCA pilot.
qZlp = reshape(out.B_ZLPCA,[],1);
qZlp = qZlp/max(norm(qZlp),eps);

[Bpos,ipos] = ba_fgls_with_pilot(out.R0,out.Rlags,qPos,cfg);
[Bzlp,izlp] = ba_fgls_with_pilot(out.R0,out.Rlags,qZlp,cfg);

cand.Bpos = Bpos;
cand.Bzlp = Bzlp;
cand.ipos = ipos;
cand.izlp = izlp;
end
