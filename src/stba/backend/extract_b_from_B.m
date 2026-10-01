function b = extract_b_from_B(B)
%EXTRACT_B_FROM_B Principal spatial vector from a rank-one PSD matrix.
%
% Global complex phase is not identifiable and is immaterial because the
% subsequent variable-projection fit re-estimates complex amplitudes.

H = (B+B')/2;
[V,D] = eig(H);
[~,idx] = max(real(diag(D)));
b = V(:,idx);

% Deterministic phase convention for reproducibility.
[~,k] = max(abs(b));
if abs(b(k)) > 0
    b = b * exp(-1j*angle(b(k)));
end
b = b / norm(b);
end
