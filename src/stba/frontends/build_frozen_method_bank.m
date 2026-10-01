function bank = build_frozen_method_bank(X, cfg)
%BUILD_FROZEN_METHOD_BANK Frozen internal method bank for paper-scale tests.
%
% IMPORTANT:
%   This helper does NOT modify any estimator developed in V0-V2A.
%   It simply collects the already-frozen front ends so that every method
%   is evaluated with the SAME continuous K-source VP-NLS backend.
%
% Methods:
%   ZL-PCA : zero-lag PCA
%   SL-CFR : single nonzero-lag common-factor recovery
%   PA-MLA : phase-aligned multi-lag averaging
%   ML-CFR : unweighted multi-lag SVD common factor
%   FGLS   : one-step feasible GLS over positive lags
%   BA-POS : bias-annihilated augmented fusion, POS pilot
%   BA-ZLP : bias-annihilated augmented fusion, zero-lag PCA pilot
%   ST-BA  : frozen split-half pilot-stability selector

out = estimate_frontends(X,cfg);
[Bfg,ifg] = fgls_cfr(out.Rlags,cfg);
[Bst,ist] = stability_adaptive_ba_fgls(X,cfg);

bank.names = {'ZL-PCA','SL-CFR','PA-MLA','ML-CFR', ...
              'FGLS','BA-POS','BA-ZLP','ST-BA'};

bank.B = {out.B_ZLPCA,out.B_SL,out.B_PA,out.B_ML, ...
          Bfg,ist.Bpos,ist.Bzlp,Bst};

bank.frontends = out;
bank.fgls = ifg;
bank.stability = ist;
end
