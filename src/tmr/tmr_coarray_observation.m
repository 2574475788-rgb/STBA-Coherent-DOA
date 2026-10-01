function out = tmr_coarray_observation(Rdec,p,duplicateLagMode)
%TMR_COARRAY_OBSERVATION Build Eqs. (20)-(23) inputs.
%
% Difference convention:
%   lag = p_i - p_j
% because R(i,j) contains exp(-j*pi*(p_i-p_j)*sin(theta)).
%
% The paper mentions redundant-data removal/correlation sorting but does
% not specify the duplicate-lag collapse rule. Default: arithmetic mean of
% all R(i,j) sharing a lag.

if nargin < 3 || isempty(duplicateLagMode)
    duplicateLagMode = 'mean';
end

p = p(:).';
M = numel(p);

allLags = zeros(M*M,1);
allVals = complex(zeros(M*M,1));
q = 0;

for j=1:M
    for i=1:M
        q=q+1;
        allLags(q)=p(i)-p(j);
        allVals(q)=Rdec(i,j);
    end
end

D = unique(allLags).';
dmin = min(D);
dmax = max(D);
Vlags = dmin:dmax;
assert(dmin == -dmax,'Expected symmetric difference coarray.');

lagVals = complex(zeros(size(Vlags)));
observed = false(size(Vlags));
multiplicity = zeros(size(Vlags));

for k=1:numel(Vlags)
    d=Vlags(k);
    idx=find(allLags==d);
    if ~isempty(idx)
        observed(k)=true;
        multiplicity(k)=numel(idx);

        switch lower(duplicateLagMode)
            case 'mean'
                lagVals(k)=mean(allVals(idx));
            case 'first'
                lagVals(k)=allVals(idx(1));
            otherwise
                error('Unknown duplicateLagMode: %s',duplicateLagMode);
        end
    else
        lagVals(k)=0;
    end
end

% Enforce the Hermitian lag symmetry implied by a Hermitian covariance.
for d=1:dmax
    ip=find(Vlags== d,1);
    im=find(Vlags==-d,1);
    if observed(ip) && observed(im)
        z=0.5*(lagVals(ip)+conj(lagVals(im)));
        lagVals(ip)=z;
        lagVals(im)=conj(z);
    end
end
i0=find(Vlags==0,1);
lagVals(i0)=real(lagVals(i0));

n=(numel(Vlags)+1)/2;  % Eq. (23), equals physical aperture+1

% Eq. (23): first column is lags 0,+1,...,+P; first row 0,-1,...,-P.
c=complex(zeros(n,1));
r=complex(zeros(1,n));
for k=0:n-1
    c(k+1)=lagVals(Vlags== k);
    r(k+1)=lagVals(Vlags==-k);
end
T=toeplitz(c,r);
T=(T+T')/2;

% Binary mask: an entry is observed iff its ULA lag exists in D.
mask=zeros(n,n);
for row=1:n
    for col=1:n
        d=(row-1)-(col-1);
        mask(row,col)=ismember(d,D);
    end
end

out.D=D;
out.Vlags=Vlags;
out.lagVals=lagVals;
out.observed=observed;
out.multiplicity=multiplicity;
out.T=T;
out.mask=mask;
out.n=n;
out.holes=Vlags(~observed);
end
