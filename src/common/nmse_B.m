function value = nmse_B(Bhat, Btrue)
%NMSE_B Normalized Frobenius error for trace-normalized spatial factors.

Bhat = Bhat / real(trace(Bhat));
Btrue = Btrue / real(trace(Btrue));
value = norm(Bhat-Btrue,'fro')^2 / norm(Btrue,'fro')^2;
end
