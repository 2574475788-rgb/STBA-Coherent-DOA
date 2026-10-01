clear; clc; close all;
thisFile=mfilename('fullpath');
root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project();
outDir=fullfile(root,'results','frozen');
if exist(outDir,'dir')~=7,mkdir(outDir);end

cfg=config_v0();
cfg.K=2;
sc=robustness_scenarios();

% ============================================================
% V3C FORMAL MODEL-MISMATCH ROBUSTNESS
%
% PRE-REGISTERED before seeing formal results.
%
% Frozen condition:
%   K=2, N=100, SNR=-10 dB.
%
% Factorial mismatch:
%   temporal = {AR1 rho=.8, fixed 7-tap FIR}
%   spatial gamma = {0, .3, .6}
%
% Same MC seed is used across all six scenarios, yielding paired geometry,
% coherent phases, innovation pool and raw spatial-noise realization.
%
% Primary internal metrics:
%   B-NMSE, DOA RMSE, P0.5.
%
% Supporting:
%   median/P95/P99 max error, ST ZLP-selection rate.
%
% External recent baseline:
%   frozen SS-corrected TMR-VAI specialization.
%
% NO algorithm or baseline parameter may be changed after viewing results.
% ============================================================

NMC=300;
N=100;
snrDb=-10;

L=8;
beamWidth=80;
minStartDist=1.5;

methodNames={'ZL-PCA','SL-CFR','PA-MLA','ML-CFR', ...
             'FGLS','BA-POS','BA-ZLP','ST-BA'};
nMethods=numel(methodNames);
nSc=numel(sc);

nmseStore=nan(NMC,nSc,nMethods);
doaMseStore=nan(NMC,nSc,nMethods);
maxErrStore=nan(NMC,nSc,nMethods);
stSelectZlp=false(NMC,nSc);

tmrMse=nan(NMC,nSc);
tmrMax=nan(NMC,nSc);

sourceGrid=cfg.angleMin:0.5:cfg.angleMax;
musicGrid=cfg.angleMin:0.1:cfg.angleMax;
dupMode='mean';

fprintf('============================================================\n');
fprintf(' V3C FORMAL MODEL-MISMATCH ROBUSTNESS\n');
fprintf(' NMC=%d, K=2, N=%d, SNR=%+.1f dB\n',NMC,N,snrDb);
fprintf(' temporal x spatial = 2 x 3 pre-registered factorial design\n');
fprintf('============================================================\n\n');

for is=1:nSc
    fprintf('%s\n',sc(is).name);

    for mc=1:NMC
        seed=cfg.seed+3320000+mc;

        [X,truth]=generate_trial_mismatch_k( ...
            cfg,N,snrDb,sc(is),[],seed);

        bank=build_frozen_method_bank(X,cfg);
        met=evaluate_method_bank_k( ...
            bank,truth,cfg,2,L,beamWidth,minStartDist);

        nmseStore(mc,is,:)=met.nmse;
        doaMseStore(mc,is,:)=met.doaMse;
        maxErrStore(mc,is,:)=met.maxErr;
        stSelectZlp(mc,is)=strcmp(bank.stability.selected,'ZLP');

        [thT,~]=tmr_vai_fully_augmented_estimate( ...
            X,cfg.p,2,sourceGrid,musicGrid,dupMode);
        [rT,mxT]=angle_error_k(thT,truth.thetaDeg);
        tmrMse(mc,is)=rT^2;
        tmrMax(mc,is)=mxT;

        if mod(mc,50)==0
            fprintf('  completed %d/%d\n',mc,NMC);
        end
    end

    checkpoint.completedIndex=is;
    checkpoint.cfg=cfg;
    checkpoint.scenarios=sc;
    checkpoint.NMC=NMC;
    checkpoint.N=N;
    checkpoint.snrDb=snrDb;
    checkpoint.methodNames=methodNames;
    checkpoint.nmseStore=nmseStore;
    checkpoint.doaMseStore=doaMseStore;
    checkpoint.maxErrStore=maxErrStore;
    checkpoint.stSelectZlp=stSelectZlp;
    checkpoint.tmrMse=tmrMse;
    checkpoint.tmrMax=tmrMax;
    save(fullfile(outDir,'v3c_robustness_checkpoint.mat'),'checkpoint','-v7.3');
end

meanNmseDb=squeeze(10*log10(mean(nmseStore,1)));
doaRmse=squeeze(sqrt(mean(doaMseStore,1)));
successProb=squeeze(mean(maxErrStore<=0.5,1));
medianMax=squeeze(median(maxErrStore,1));
p95Max=squeeze(prctile(maxErrStore,95,1));
p99Max=squeeze(prctile(maxErrStore,99,1));
stRate=mean(stSelectZlp,1);

tmrRmse=sqrt(mean(tmrMse,1));
tmrP=mean(tmrMax<=0.5,1);
tmrMed=median(tmrMax,1);
tmrP95=prctile(tmrMax,95,1);
tmrP99=prctile(tmrMax,99,1);

idxFG=5; idxPOS=6; idxZLP=7; idxST=8;
Bboot=3000;

fprintf('\n==================== V3C FORMAL RESULTS ====================\n');

for is=1:nSc
    fprintf('\n%-10s | ST ZLP-select=%.3f\n',sc(is).name,stRate(is));

    for m=[idxFG idxPOS idxZLP idxST]
        fprintf('  %-7s B-NMSE=%8.3f dB | RMSE=%8.4f | P0.5=%6.3f | med=%7.4f | P95=%7.4f | P99=%7.4f\n', ...
            methodNames{m},meanNmseDb(is,m),doaRmse(is,m), ...
            successProb(is,m),medianMax(is,m),p95Max(is,m),p99Max(is,m));
    end

    fprintf('  TMR-VAI          RMSE=%8.4f | P0.5=%6.3f | med=%7.4f | P95=%7.4f | P99=%7.4f\n', ...
        tmrRmse(is),tmrP(is),tmrMed(is),tmrP95(is),tmrP99(is));

    % Internal paired comparisons
    comps={ ...
        'ST-FG',idxST,idxFG; ...
        'ST-POS',idxST,idxPOS; ...
        'ST-ZLP',idxST,idxZLP};

    for cc=1:size(comps,1)
        label=comps{cc,1}; ia=comps{cc,2}; ib=comps{cc,3};

        dN=nmseStore(:,is,ia)-nmseStore(:,is,ib);
        dD=doaMseStore(:,is,ia)-doaMseStore(:,is,ib);
        dP=double(maxErrStore(:,is,ia)<=0.5) - ...
           double(maxErrStore(:,is,ib)<=0.5);

        [loN,hiN]=paired_bootstrap_ci(dN,Bboot,20310000+100*is+cc);
        [loD,hiD]=paired_bootstrap_ci(dD,Bboot,20311000+100*is+cc);
        [loP,hiP]=paired_bootstrap_ci(dP,Bboot,20312000+100*is+cc);

        fprintf('  %-6s dNMSE(linear)=%+10.4e CI[%+10.4e,%+10.4e] | dDOA-MSE=%+10.4e CI[%+10.4e,%+10.4e] | dP0.5=%+7.4f CI[%+7.4f,%+7.4f]\n', ...
            label,mean(dN),loN,hiN,mean(dD),loD,hiD,mean(dP),loP,hiP);
    end

    % External paired ST vs TMR
    dDT=doaMseStore(:,is,idxST)-tmrMse(:,is);
    dPT=double(maxErrStore(:,is,idxST)<=0.5) - ...
        double(tmrMax(:,is)<=0.5);

    [loDT,hiDT]=paired_bootstrap_ci(dDT,Bboot,20313000+is);
    [loPT,hiPT]=paired_bootstrap_ci(dPT,Bboot,20314000+is);

    fprintf('  ST-TMR dDOA-MSE=%+10.4e CI[%+10.4e,%+10.4e] | dP0.5=%+7.4f CI[%+7.4f,%+7.4f]\n', ...
        mean(dDT),loDT,hiDT,mean(dPT),loPT,hiPT);
end

%% Paired mismatch penalties for ST relative to in-model AR1-WHITE
fprintf('\n==================== ST-BA PAIRED MISMATCH PENALTIES ====================\n');

for is=2:nSc
    dN=nmseStore(:,is,idxST)-nmseStore(:,1,idxST);
    dD=doaMseStore(:,is,idxST)-doaMseStore(:,1,idxST);
    dP=double(maxErrStore(:,is,idxST)<=0.5) - ...
       double(maxErrStore(:,1,idxST)<=0.5);

    [loN,hiN]=paired_bootstrap_ci(dN,Bboot,20315000+is);
    [loD,hiD]=paired_bootstrap_ci(dD,Bboot,20316000+is);
    [loP,hiP]=paired_bootstrap_ci(dP,Bboot,20317000+is);

    fprintf('%-10s vs AR1-WHITE: dNMSE=%+10.4e CI[%+10.4e,%+10.4e] | dDOA-MSE=%+10.4e CI[%+10.4e,%+10.4e] | dP0.5=%+7.4f CI[%+7.4f,%+7.4f]\n', ...
        sc(is).name,mean(dN),loN,hiN,mean(dD),loD,hiD,mean(dP),loP,hiP);
end

results.cfg=cfg;
results.scenarios=sc;
results.NMC=NMC;
results.N=N;
results.snrDb=snrDb;
results.methodNames=methodNames;
results.nmseStore=nmseStore;
results.doaMseStore=doaMseStore;
results.maxErrStore=maxErrStore;
results.stSelectZlp=stSelectZlp;
results.meanNmseDb=meanNmseDb;
results.doaRmse=doaRmse;
results.successProb=successProb;
results.medianMax=medianMax;
results.p95Max=p95Max;
results.p99Max=p99Max;
results.stRate=stRate;
results.tmrMse=tmrMse;
results.tmrMax=tmrMax;
results.tmrRmse=tmrRmse;
results.tmrP=tmrP;
results.tmrMed=tmrMed;
results.tmrP95=tmrP95;
results.tmrP99=tmrP99;

save(fullfile(outDir,'v3c_robustness_formal_results.mat'),'results','-v7.3');

figure;
plot(1:nSc,doaRmse(:,idxST),'-o','LineWidth',1.1); hold on;
plot(1:nSc,tmrRmse,'-s','LineWidth',1.1);
grid on; xticks(1:nSc); xticklabels({sc.name}); xtickangle(30);
ylabel('DOA RMSE (deg)'); legend('ST-BA','TMR-VAI','Location','best');
title('Model-mismatch robustness: DOA RMSE');

figure;
plot(1:nSc,successProb(:,idxST),'-o','LineWidth',1.1); hold on;
plot(1:nSc,tmrP,'-s','LineWidth',1.1);
grid on; xticks(1:nSc); xticklabels({sc.name}); xtickangle(30);
ylabel('P(max error \leq 0.5^\circ)'); ylim([0 1.02]);
legend('ST-BA','TMR-VAI','Location','best');
title('Model-mismatch robustness: success probability');

fprintf('\nSaved: %s\n',fullfile(outDir,'v3c_robustness_formal_results.mat'));
