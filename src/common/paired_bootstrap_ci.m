function [lo,hi] = paired_bootstrap_ci(x,B,seed)
%PAIRED_BOOTSTRAP_CI Percentile bootstrap CI for a paired mean difference.
x = x(:);
x = x(isfinite(x));
n = numel(x);
if n < 2
    lo = NaN; hi = NaN; return;
end

rng(seed,'twister');
bm = zeros(B,1);
for bb = 1:B
    idx = randi(n,n,1);
    bm(bb) = mean(x(idx));
end
bm = sort(bm);

lo = qsorted(bm,2.5);
hi = qsorted(bm,97.5);
end

function q = qsorted(xs,p)
n = numel(xs);
pos = 1 + (p/100)*(n-1);
i = floor(pos);
f = pos-i;
if i >= n
    q = xs(end);
else
    q = xs(i) + f*(xs(i+1)-xs(i));
end
end
