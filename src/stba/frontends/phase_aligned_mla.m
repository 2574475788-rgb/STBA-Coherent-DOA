function [B, info] = phase_aligned_mla(Rlags, cfg)
%PHASE_ALIGNED_MLA Strong multi-lag averaging baseline.
%
% Different lags have different complex correlation phases. A direct
% arithmetic average can cancel useful signal. Therefore each lag is
% phase-aligned to the first lag before averaging.
%
% c_j = <R_1,R_j>_F = trace(R_1^H R_j)
% gamma_j = c_j/|c_j|
% R_j_aligned = conj(gamma_j)*R_j

[~,~,J] = size(Rlags);
Rref = Rlags(:,:,1);
Rsum = complex(zeros(size(Rref)));
gammas = ones(J,1);

for jj = 1:J
    Rj = Rlags(:,:,jj);
    c = sum(conj(Rref(:)).*Rj(:));
    if abs(c) > cfg.epsPhase
        gamma = c/abs(c);
    else
        gamma = 1;
    end
    gammas(jj) = gamma;
    Rsum = Rsum + conj(gamma)*Rj;
end

Ravg = Rsum/J;
[B, pinfo] = project_rank1_psd(Ravg, cfg);

info.gammas = gammas;
info.Ravg = Ravg;
info.projection = pinfo;
end
