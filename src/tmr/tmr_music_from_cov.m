function [thetaHat,P,info] = tmr_music_from_cov(R,K,thetaGrid)
%TMR_MUSIC_FROM_COV Known-K MUSIC on reconstructed ULA covariance.

R=(R+R')/2;
n=size(R,1);

[U,D]=eig(R,'vector');
[D,ord]=sort(real(D),'descend');
U=U(:,ord);

if K >= n
    error('MUSIC requires K < covariance dimension.');
end

Un=U(:,K+1:end);
Pn=Un*Un';

A=tmr_steering(0:n-1,thetaGrid);
den=real(sum(conj(A).*(Pn*A),1));
den=max(den,eps);
P=1./den;

thetaHat=pick_peaks(thetaGrid,P,K);

info.eigenvalues=D;
info.noiseSubspace=Un;
end

function thetaHat=pick_peaks(grid,P,K)
grid=grid(:).';
P=P(:).';

isPeak=false(size(P));
if numel(P)>=3
    isPeak(2:end-1)=P(2:end-1)>=P(1:end-2) & P(2:end-1)>=P(3:end);
end
isPeak(1)=P(1)>=P(2);
isPeak(end)=P(end)>=P(end-1);

cand=find(isPeak);
[~,ord]=sort(P(cand),'descend');
cand=cand(ord);

sel=[];
for ii=1:numel(cand)
    sel(end+1)=cand(ii); %#ok<AGROW>
    if numel(sel)==K, break; end
end

% Defensive fallback if a very flat spectrum yields too few local maxima.
if numel(sel)<K
    [~,ordAll]=sort(P,'descend');
    for ii=1:numel(ordAll)
        if ~ismember(ordAll(ii),sel)
            sel(end+1)=ordAll(ii); %#ok<AGROW>
            if numel(sel)==K, break; end
        end
    end
end

thetaHat=sort(grid(sel));
end
