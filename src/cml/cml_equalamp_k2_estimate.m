function [thetaHat, info] = cml_equalamp_k2_estimate(X, p, opts)
%CML_EQUALAMP_K2_ESTIMATE Source-faithful equal-amplitude K=2 CML baseline.
%
% Implements the concentrated deterministic maximum-likelihood benchmark
% described in S. Zhang et al., Digital Signal Processing 184 (2026) 106468,
% Sec. 3.4 / simulation implementation:
%
%   b(theta1,theta2,phi) = a(theta1) + exp(1j*phi) a(theta2)
%   J = (b^H Rhat b)/(b^H b),  theta1 < theta2.
%
% Default search protocol follows that paper exactly:
%   coarse DOA step       0.5 deg over [-60,60]
%   coarse phase step    10   deg over [-180,180)
%   retain                5 spatially distinct coarse solutions
%   coarse distinctness   0.75 deg Euclidean distance in [theta1 theta2]
%   fine DOA window      +/-0.5 deg, step 0.1 deg
%   fine phase window    +/-10 deg, step 1 deg
%
% The likelihood formula itself is valid for an arbitrary known linear-array
% sensor-position vector p (in half-wavelength units); no ULA structure is
% used by this implementation. This makes the benchmark directly applicable
% to the fixed nested sparse array used by ST-BA without altering the CML
% statistical model.
%
% IMPORTANT SCOPE:
%   * K is fixed to 2.
%   * Coherent coefficients have equal magnitude [1, exp(j*phi)].
%   * The true phase is NOT an input and is never used by the estimator.
%   * This is an accuracy-oriented likelihood benchmark, not a low-complexity
%     method.

if nargin < 3 || isempty(opts), opts = struct(); end
opts = fill_defaults(opts);

validateattributes(X,{'numeric'},{'2d','nonempty','finite'});
validateattributes(p,{'numeric'},{'vector','real','finite'});
assert(size(X,1)==numel(p),'Rows of X must match numel(p).');
assert(opts.angleMin < opts.angleMax,'angleMin must be < angleMax.');
assert(opts.coarseAngleStep > 0 && opts.fineAngleStep > 0,'Angle steps must be positive.');
assert(opts.coarsePhaseStep > 0 && opts.finePhaseStep > 0,'Phase steps must be positive.');
assert(opts.retainCount >= 1 && opts.retainCount==floor(opts.retainCount),'retainCount must be a positive integer.');

N = size(X,2);
Rx = (X*X')/N;
Rx = (Rx+Rx')/2;
M = numel(p);

%% Coarse global search
coarseAngles = make_grid(opts.angleMin,opts.coarseAngleStep,opts.angleMax);
coarsePhases = -180:opts.coarsePhaseStep:(180-opts.coarsePhaseStep);
A = steering_matrix(p,coarseAngles);
Q = A'*Rx*A;
G = A'*A;

nG = numel(coarseAngles);
[ii,jj] = find(triu(true(nG),1));
qii = real(Q(sub2ind([nG nG],ii,ii)));
qjj = real(Q(sub2ind([nG nG],jj,jj)));
qij = Q(sub2ind([nG nG],ii,jj));
gij = G(sub2ind([nG nG],ii,jj));

C = exp(1j*deg2rad(coarsePhases));
num = qii + qjj + 2*real(qij .* C);
den = 2*M + 2*real(gij .* C);
score = num ./ den;
score(~isfinite(score) | den <= opts.denomTol*M) = -Inf;

[sortedScore,linearOrder] = sort(score(:),'descend');
retained = zeros(opts.retainCount,4); % theta1 theta2 phi score
nRet = 0;
for kk=1:numel(linearOrder)
    if ~isfinite(sortedScore(kk)), break; end
    [pairIdx,phaseIdx] = ind2sub(size(score),linearOrder(kk));
    candAngles = [coarseAngles(ii(pairIdx)), coarseAngles(jj(pairIdx))];
    if nRet > 0
        d = sqrt(sum((retained(1:nRet,1:2)-candAngles).^2,2));
        if any(d < opts.coarseDistinctDist), continue; end
    end
    nRet = nRet+1;
    retained(nRet,:) = [candAngles, coarsePhases(phaseIdx), sortedScore(kk)];
    if nRet == opts.retainCount, break; end
end
assert(nRet>=1,'CML coarse search produced no finite candidate.');
retained = retained(1:nRet,:);

%% Local fine searches around retained coarse solutions
bestScore = -Inf;
best = [NaN NaN NaN];
fineHypCount = 0;
for rr=1:nRet
    t1c = retained(rr,1); t2c = retained(rr,2); pc = retained(rr,3);
    t1Grid = make_grid(max(opts.angleMin,t1c-opts.fineAngleHalfWidth), ...
        opts.fineAngleStep,min(opts.angleMax,t1c+opts.fineAngleHalfWidth));
    t2Grid = make_grid(max(opts.angleMin,t2c-opts.fineAngleHalfWidth), ...
        opts.fineAngleStep,min(opts.angleMax,t2c+opts.fineAngleHalfWidth));
    phiGrid = pc + (-opts.finePhaseHalfWidth:opts.finePhaseStep:opts.finePhaseHalfWidth);
    phiGrid = wrap_phase_deg(phiGrid);

    [T1,T2] = ndgrid(t1Grid,t2Grid);
    valid = T1 < T2;
    t1 = T1(valid); t2 = T2(valid);
    if isempty(t1), continue; end

    A1 = steering_matrix(p,t1(:).');
    A2 = steering_matrix(p,t2(:).');
    q11 = real(sum(conj(A1).*(Rx*A1),1)).';
    q22 = real(sum(conj(A2).*(Rx*A2),1)).';
    q12 = sum(conj(A1).*(Rx*A2),1).';
    g12 = sum(conj(A1).*A2,1).';

    Cf = exp(1j*deg2rad(phiGrid));
    numf = q11 + q22 + 2*real(q12 .* Cf);
    denf = 2*M + 2*real(g12 .* Cf);
    sf = numf ./ denf;
    sf(~isfinite(sf) | denf <= opts.denomTol*M) = -Inf;
    fineHypCount = fineHypCount + numel(sf);

    [smax,idx] = max(sf(:));
    if isfinite(smax) && smax > bestScore
        [pairIdx,phaseIdx] = ind2sub(size(sf),idx);
        bestScore = smax;
        best = [t1(pairIdx), t2(pairIdx), phiGrid(phaseIdx)];
    end
end

assert(all(isfinite(best)),'CML fine search produced no finite candidate.');
thetaHat = sort(best(1:2));

% The fine search includes every retained coarse center. Thus its optimum
% should never be worse than the best retained coarse score, up to roundoff.
bestCoarseScore = max(retained(:,4));
if bestScore + opts.scoreTol < bestCoarseScore
    error('CML fidelity check failed: fine optimum < retained coarse optimum.');
end

info.Rx = Rx;
info.thetaHat = thetaHat;
info.phiHatDeg = best(3);
info.finalScore = bestScore;
info.bestCoarseScore = bestCoarseScore;
info.retainedCoarse = retained;
info.options = opts;
info.coarseAngleCount = nG;
info.coarsePairCount = numel(ii);
info.coarsePhaseCount = numel(coarsePhases);
info.coarseHypothesisCount = numel(score);
info.fineHypothesisCount = fineHypCount;
info.source = 'S. Zhang et al., Digital Signal Processing 184 (2026) 106468, CML benchmark';
info.truePhaseUsed = false;
end

function opts = fill_defaults(opts)
def.angleMin = -60;
def.angleMax = 60;
def.coarseAngleStep = 0.5;
def.coarsePhaseStep = 10;
def.retainCount = 5;
def.coarseDistinctDist = 0.75;
def.fineAngleHalfWidth = 0.5;
def.fineAngleStep = 0.1;
def.finePhaseHalfWidth = 10;
def.finePhaseStep = 1;
def.denomTol = 1e-12;
def.scoreTol = 1e-9;
fn = fieldnames(def);
for k=1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = def.(fn{k});
    end
end
end

function g = make_grid(lo,step,hi)
% Stable inclusive decimal grid without accumulated endpoint drift.
n = floor((hi-lo)/step + 1e-10);
g = lo + (0:n)*step;
if isempty(g) || hi-g(end) > 0.25*step
    g(end+1) = hi;
elseif abs(g(end)-hi) <= 1e-9
    g(end)=hi;
end
end

function p = wrap_phase_deg(p)
p = mod(p+180,360)-180;
% Keep stable sorted unique values; -180 and +180 represent same phase.
p = unique(round(p*1e12)/1e12,'stable');
end
