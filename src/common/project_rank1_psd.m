function [B, info] = project_rank1_psd(Q, cfg)
%PROJECT_RANK1_PSD Recover a Hermitian PSD rank-one spatial factor.
%
% Steps:
%   1) fix the arbitrary global complex phase using trace(Q);
%   2) Hermitian projection;
%   3) retain the largest eigenvector (rank-one PSD projection);
%   4) trace-normalize.
%
% The scale is immaterial for the subsequent DOA diagnostic.

trQ = trace(Q);

if abs(trQ) > cfg.epsPhase
    Q0 = exp(-1j*angle(trQ)) * Q;
else
    % This should not occur for a meaningful common factor.
    Q0 = Q;
end

H = (Q0 + Q0')/2;
[V,D] = eig(H);
evals = real(diag(D));
[lambda1, idx] = max(evals);
u1 = V(:,idx);

if lambda1 <= cfg.epsEig
    error('project_rank1_psd:NonPositive', ...
        'Largest Hermitian eigenvalue is non-positive (%.3e).', lambda1);
end

% Since final scale is normalized, use u1*u1^H directly.
B = u1*u1';
B = (B + B')/2;
B = B / real(trace(B));

info.traceBefore = trQ;
info.lambda1 = lambda1;
info.eigvals = evals;
info.hermitianError = norm(H-H','fro') / max(norm(H,'fro'),eps);
end
