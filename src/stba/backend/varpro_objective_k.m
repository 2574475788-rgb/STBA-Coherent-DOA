function f = varpro_objective_k(theta, b, cfg, K)
%VARPRO_OBJECTIVE_K K-source separable nonlinear least-squares objective.
%
% From the recovered coherent spatial factor:
%   b ~= A(theta) beta,
%
% eliminate the unknown complex coherent coefficients:
%   beta_hat(theta) = argmin_beta ||b-A(theta)beta||_2.
%
% The profiled objective is
%   f(theta) = ||P_A(theta)^perp b||_2^2 / ||b||_2^2.
%
% This is the natural K-source generalization of the V0 K=2 diagnostic
% variable-projection decoder.

theta = sort(theta(:).');

if numel(theta) ~= K
    f = 1e3 + abs(numel(theta)-K);
    return;
end

% FOV guard.
if any(theta < cfg.angleMin) || any(theta > cfg.angleMax)
    low = max(cfg.angleMin-theta,0);
    high = max(theta-cfg.angleMax,0);
    f = 1e3 + sum(low.^2 + high.^2);
    return;
end

% Near-collision guard.
if K > 1 && min(diff(theta)) < cfg.refineMinSep
    gap = min(diff(theta));
    f = 1e2 + (cfg.refineMinSep-gap)^2;
    return;
end

A = steering_matrix(cfg.p,theta);

% With the separation guard and the present nested array, A is normally
% full column rank. Backslash gives the exact LS variable-projection step.
betaHat = A\b;
r = b-A*betaHat;

f = real(r'*r)/max(real(b'*b),eps);

if ~isfinite(f)
    f = 1e6;
end
end
