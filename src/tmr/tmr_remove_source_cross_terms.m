function Rdec = tmr_remove_source_cross_terms(Rx,A,Rs)
%TMR_REMOVE_SOURCE_CROSS_TERMS Operational form of Eqs. (18)-(19).
%
% The paper states that the off-diagonal source-covariance contribution
% A * (Rs-diag(diag(Rs))) * A^H is subtracted.
%
% Therefore, using the sample covariance directly:
%
%   Rdec = Rx - A*offdiag(Rs)*A^H.
%
% This avoids inventing a separate noise-variance estimator, which is not
% specified in Algorithm 1.

off = Rs - diag(diag(Rs));
Rdec = Rx - A*off*A';
Rdec = (Rdec+Rdec')/2;
end
