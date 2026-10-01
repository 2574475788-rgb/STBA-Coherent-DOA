clear; clc; close all;
thisFile=mfilename('fullpath');
root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project();
outDir=fullfile(root,'results','frozen');
if exist(outDir,'dir')~=7,mkdir(outDir);end
cfg = config_v0();

% ============================================================
% Version 1G: split-half pilot-stability selector
%
% V1F result:
%   cross-fitting corrected the strong high-rho ZLP over-selection:
%       rho=1.00  IN=0.60 -> CF=0.27
%       rho=0.95  IN=0.83 -> CF=0.37
%   and preserved the low-rho transition.
%
% But per-trial oracle agreement remained only moderate, and CF can still
% choose the wrong branch in low-rho trials where the cost is large.
%
% V1G therefore stops scoring FINAL BA candidates and instead directly
% measures the quantity V1D identified as the bottleneck: PILOT stability.
%
% Selector:
%   kappa_POS = split-half Frobenius consistency of POS pilots
%   kappa_ZLP = split-half Frobenius consistency of ZLP pilots
%   choose the larger kappa.
%
% Front-end only. DOA decoding is intentionally omitted.
% ============================================================

NMC = 100;
N = 100;
snrDb = -10;
rhoList = [1.00 0.95 0.90 0.80 0.70];

methodNames = {'ZL-PCA','ML-CFR','BA-POS','BA-ZLP','CF-BA','ST-BA'};
nMethods = numel(methodNames);
nRho = numel(rhoList);

nmseStore = nan(NMC,nRho,nMethods);

cfSelectZlp = false(NMC,nRho);
stSelectZlp = false(NMC,nRho);
oracleZlp = false(NMC,nRho);

kPosStore = nan(NMC,nRho);
kZlpStore = nan(NMC,nRho);

oracleSelectorNmse = nan(NMC,nRho);

fprintf('============================================================\n');
fprintf(' V1G split-half pilot-stability selector\n');
fprintf(' NMC=%d, N=%d, SNR=%g dB\n',NMC,N,snrDb);
fprintf(' Select larger split-half pilot consistency; no threshold/rho\n');
fprintf('============================================================\n\n');

for mc=1:NMC
    seed = cfg.seed + 1400000 + mc;

    for ir=1:nRho
        rho = rhoList(ir);

        [X,truth] = generate_trial_ar1(cfg,N,snrDb,rho,[],seed);

        out = estimate_frontends(X,cfg);

        % Previous cross-fitted selector.
        [Bcf,icf] = crossfit_adaptive_ba_fgls(X,cfg);

        % New stability selector.
        [Bst,ist] = stability_adaptive_ba_fgls(X,cfg);

        % Use the same full-data branch candidates returned by ST selector.
        Bpos = ist.Bpos;
        Bzlp = ist.Bzlp;

        Blist = {out.B_ZLPCA,out.B_ML,Bpos,Bzlp,Bcf,Bst};

        for m=1:nMethods
            nmseStore(mc,ir,m) = nmse_B(Blist{m},truth.Btrue);
        end

        cfSelectZlp(mc,ir) = strcmp(icf.selected,'ZLP');
        stSelectZlp(mc,ir) = strcmp(ist.selected,'ZLP');

        kPosStore(mc,ir) = ist.kappaPos;
        kZlpStore(mc,ir) = ist.kappaZlp;

        nPos = nmse_B(Bpos,truth.Btrue);
        nZlp = nmse_B(Bzlp,truth.Btrue);
        oracleZlp(mc,ir) = nZlp<nPos;
        oracleSelectorNmse(mc,ir) = min(nPos,nZlp);
    end

    if mod(mc,10)==0
        fprintf('Completed %d/%d trials\n',mc,NMC);
    end
end

meanNmseDb = squeeze(10*log10(mean(nmseStore,1)));

cfRate = mean(cfSelectZlp,1);
stRate = mean(stSelectZlp,1);
orRate = mean(oracleZlp,1);

cfAgree = mean(cfSelectZlp==oracleZlp,1);
stAgree = mean(stSelectZlp==oracleZlp,1);

meanKPos = mean(kPosStore,1);
meanKZlp = mean(kZlpStore,1);

oracleSelectorDb = 10*log10(mean(oracleSelectorNmse,1));

idxZL=1; idxML=2; idxPOS=3; idxZLP=4; idxCF=5; idxST=6;
Bboot=2000;

fprintf('\n==================== V1G RESULTS ====================\n');
for ir=1:nRho
    rho = rhoList(ir);

    fprintf('\nrho = %.2f\n',rho);
    fprintf('  mean pilot stability: POS=%.4f, ZLP=%.4f, margin(ZLP-POS)=%+.4f\n', ...
        meanKPos(ir),meanKZlp(ir),meanKZlp(ir)-meanKPos(ir));
    fprintf('  ZLP selection rate  : CF=%.3f, ST=%.3f, ORACLE=%.3f\n', ...
        cfRate(ir),stRate(ir),orRate(ir));
    fprintf('  oracle agreement    : CF=%.3f, ST=%.3f\n', ...
        cfAgree(ir),stAgree(ir));

    for m=1:nMethods
        fprintf('  %-7s B-NMSE=%8.3f dB\n',methodNames{m},meanNmseDb(ir,m));
    end
    fprintf('  OR-SEL  B-NMSE=%8.3f dB [diagnostic best-of-two]\n',oracleSelectorDb(ir));

    comps = { ...
        'ST-CF',idxST,idxCF; ...
        'ST-POS',idxST,idxPOS; ...
        'ST-ZLP',idxST,idxZLP; ...
        'ST-ML',idxST,idxML; ...
        'ST-ZL',idxST,idxZL};

    for cc=1:size(comps,1)
        label=comps{cc,1};
        ia=comps{cc,2};
        ib=comps{cc,3};

        d=nmseStore(:,ir,ia)-nmseStore(:,ir,ib);
        [lo,hi]=paired_bootstrap_ci(d,Bboot,20272000+100*ir+cc);

        fprintf('  %-6s delta linear B-NMSE=%+10.4e CI[%+10.4e,%+10.4e] win=%.3f\n', ...
            label,mean(d),lo,hi,mean(d<0));
    end
end

results.cfg=cfg;
results.NMC=NMC;
results.N=N;
results.snrDb=snrDb;
results.rhoList=rhoList;
results.methodNames=methodNames;
results.nmseStore=nmseStore;
results.meanNmseDb=meanNmseDb;
results.cfRate=cfRate;
results.stRate=stRate;
results.oracleRate=orRate;
results.cfAgreement=cfAgree;
results.stAgreement=stAgree;
results.meanKPos=meanKPos;
results.meanKZlp=meanKZlp;
results.oracleSelectorDb=oracleSelectorDb;

save(fullfile(outDir,'v1g_stability_selector_results.mat'),'results','-v7.3');

figure;
plot(rhoList,[meanKPos(:),meanKZlp(:)],'-o','LineWidth',1.2);
grid on;
xlabel('\rho');
ylabel('Mean split-half pilot consistency');
ylim([0 1.02]);
legend({'POS pilot','ZLP pilot'},'Location','best');
title('V1G pilot stability');
set(gca,'XDir','reverse');

figure;
plot(rhoList,[cfRate(:),stRate(:),orRate(:)],'-o','LineWidth',1.2);
grid on;
xlabel('\rho');
ylabel('ZLP selection rate');
ylim([0 1.05]);
legend({'Cross-fit score','Stability','Oracle'},'Location','best');
title('V1G selector behavior');
set(gca,'XDir','reverse');

figure;
plot(rhoList,meanNmseDb,'-o','LineWidth',1.2);
grid on;
xlabel('\rho');
ylabel('Mean B-NMSE (dB)');
legend(methodNames,'Location','best');
title(sprintf('V1G front-end recovery, SNR=%g dB, N=%d',snrDb,N));
set(gca,'XDir','reverse');

fprintf('\nSaved: %s\n',fullfile(outDir,'v1g_stability_selector_results.mat'));
fprintf('Decision focus: does pilot stability improve per-trial oracle agreement and low-rho branch selection?\n');
