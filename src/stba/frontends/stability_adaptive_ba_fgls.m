function [B, info] = stability_adaptive_ba_fgls(X, cfg)
%STABILITY_ADAPTIVE_BA_FGLS Split-half pilot-stability selector.
%
% V1D showed that the main bottleneck is PILOT quality:
%   high temporal correlation -> POS pilot is stronger;
%   low temporal correlation  -> ZLP pilot is stronger.
%
% Instead of scoring the final BA candidates with the same likelihood model,
% this selector directly asks which PILOT family is more reproducible across
% two non-overlapping temporal halves.
%
% Split:
%   A = first half, B = second half.
%
% For each family f in {POS,ZLP}, compute structured pilot matrices
%   B_f^A, B_f^B
% and their normalized Frobenius consistency
%
%   kappa_f =
%     |<B_f^A,B_f^B>_F| /
%     (||B_f^A||_F ||B_f^B||_F),      0 <= kappa_f <= 1.
%
% Select the family with the larger kappa_f.
% Then recompute that family on ALL snapshots and apply the SAME BA-FGLS
% update used previously.
%
% No rho, ground truth, empirical threshold, or tuning parameter is used.

N = size(X,2);
nA = floor(N/2);

assert(nA > max(cfg.lags) && (N-nA) > max(cfg.lags), ...
    'Each half must contain more samples than max(cfg.lags).');

XA = X(:,1:nA);
XB = X(:,nA+1:end);

% Half-data raw pilots.
[pA,outA] = make_raw_pilots(XA,cfg);
[pB,outB] = make_raw_pilots(XB,cfg);

kPos = fro_consistency(pA.Bpos,pB.Bpos);
kZlp = fro_consistency(pA.Bzlp,pB.Bzlp);

% Full-data pilots / observations.
[pF,outF] = make_raw_pilots(X,cfg);

% Compute both final candidates for diagnostics and consistent comparison.
[Bpos,ipos] = ba_fgls_with_pilot(outF.R0,outF.Rlags,pF.qPos,cfg);
[Bzlp,izlp] = ba_fgls_with_pilot(outF.R0,outF.Rlags,pF.qZlp,cfg);

if kZlp > kPos
    B = Bzlp;
    selected = 'ZLP';
    selectedIndex = 2;
else
    B = Bpos;
    selected = 'POS';
    selectedIndex = 1;
end

info.selected = selected;
info.selectedIndex = selectedIndex;
info.kappaPos = kPos;
info.kappaZlp = kZlp;
info.kappaMargin = kZlp-kPos; % positive => ZLP
info.Bpos = Bpos;
info.Bzlp = Bzlp;
info.pos = ipos;
info.zlp = izlp;
info.halfA = pA;
info.halfB = pB;
info.full = pF;
info.outA = outA;
info.outB = outB;
info.outFull = outF;
end

function [p,out] = make_raw_pilots(X,cfg)
out = estimate_frontends(X,cfg);

[M,~,J] = size(out.Rlags);
Z = complex(zeros(M*M,J));
for jj=1:J
    Z(:,jj) = reshape(out.Rlags(:,:,jj),[],1);
end

% POS pilot: dominant left singular vector, then same structural projection
% used elsewhere so consistency is measured between physically comparable B.
[U,~,~] = svd(Z,'econ');
q = U(:,1);
Q = reshape(q,M,M);
[Bpos,~] = project_rank1_psd(Q,cfg);
qPos = reshape(Bpos,[],1);
qPos = qPos/max(norm(qPos),eps);

% ZLP pilot.
Bzlp = out.B_ZLPCA;
qZlp = reshape(Bzlp,[],1);
qZlp = qZlp/max(norm(qZlp),eps);

p.Bpos = Bpos;
p.Bzlp = Bzlp;
p.qPos = qPos;
p.qZlp = qZlp;
end

function k = fro_consistency(A,B)
den = norm(A,'fro')*norm(B,'fro');
if den <= eps
    k = 0;
else
    k = abs(sum(conj(A(:)).*B(:))) / den;
    k = min(max(real(k),0),1);
end
end
