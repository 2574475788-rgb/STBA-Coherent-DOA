function [thetaCand, info] = coarse_k_beam_candidates(b, cfg, K, L, beamWidth, minStartDist)
%COARSE_K_BEAM_CANDIDATES Deterministic K-source coarse initialization.
%
% The full K-dimensional grid grows combinatorially. Instead of a single
% greedy OMP path, retain a BEAM of the best partial supports:
%
%   depth 1: evaluate all grid atoms;
%   depth 2...K: expand every retained support by every admissible grid atom,
%                deduplicate, score by variable-projection residual,
%                retain the best beamWidth supports.
%
% At the final depth, return up to L spatially distinct starts.
%
% This is only an initializer. Final DOAs are continuous because every start
% is refined by K-dimensional variable projection.

if nargin < 4 || isempty(L), L = 8; end
if nargin < 5 || isempty(beamWidth), beamWidth = 80; end
if nargin < 6 || isempty(minStartDist), minStartDist = 1.5; end

grid = cfg.coarseGrid(:).';
G = numel(grid);
Agrid = steering_matrix(cfg.p,grid);

bn2 = real(b'*b);

% ----- depth 1 -----
supports = (1:G).';
scores = inf(G,1);

for ii=1:G
    a=Agrid(:,ii);
    c=(a'*b)/(a'*a);
    r=b-a*c;
    scores(ii)=real(r'*r)/max(bn2,eps);
end

[~,ord]=sort(scores,'ascend');
keep=ord(1:min(beamWidth,numel(ord)));
supports=supports(keep,:);
scores=scores(keep);

% ----- deeper supports -----
for depth=2:K
    cand = zeros(beamWidth*G,depth);
    nCand=0;

    for rr=1:size(supports,1)
        old=supports(rr,:);
        for gg=1:G
            inds=sort([old,gg]);

            % no duplicate grid atom
            if numel(unique(inds)) < depth
                continue;
            end

            th=grid(inds);
            if min(diff(th)) < cfg.coarseMinSep
                continue;
            end

            nCand=nCand+1;
            cand(nCand,:)=inds;
        end
    end

    cand=cand(1:nCand,:);
    if isempty(cand)
        error('coarse_k_beam_candidates:NoCandidates', ...
            'No admissible coarse supports at depth %d.',depth);
    end

    % Deduplicate permutations / repeated expansions.
    cand=unique(cand,'rows','stable');

    candScores=inf(size(cand,1),1);
    for cc=1:size(cand,1)
        A=Agrid(:,cand(cc,:));
        x=A\b;
        r=b-A*x;
        candScores(cc)=real(r'*r)/max(bn2,eps);
    end

    [~,ord]=sort(candScores,'ascend');
    keep=ord(1:min(beamWidth,numel(ord)));
    supports=cand(keep,:);
    scores=candScores(keep);
end

% ----- choose spatially distinct final starts -----
[~,ord]=sort(scores,'ascend');
supports=supports(ord,:);
scores=scores(ord);

thetaCand=zeros(0,K);
scoreCand=zeros(0,1);
supportCand=zeros(0,K);

for rr=1:size(supports,1)
    th=grid(supports(rr,:));

    if isempty(thetaCand)
        accept=true;
    else
        d=sqrt(sum((thetaCand-th).^2,2));
        accept=all(d >= minStartDist);
    end

    if accept
        thetaCand(end+1,:)=th; %#ok<AGROW>
        scoreCand(end+1,1)=scores(rr); %#ok<AGROW>
        supportCand(end+1,:)=supports(rr,:); %#ok<AGROW>
    end

    if size(thetaCand,1) >= L
        break;
    end
end

% Always return at least the best support.
if isempty(thetaCand)
    thetaCand=grid(supports(1,:));
    scoreCand=scores(1);
    supportCand=supports(1,:);
end

info.grid=grid;
info.beamWidth=beamWidth;
info.L=L;
info.minStartDist=minStartDist;
info.supports=supportCand;
info.scores=scoreCand;
info.bestCoarseResidual=scoreCand(1);
end
