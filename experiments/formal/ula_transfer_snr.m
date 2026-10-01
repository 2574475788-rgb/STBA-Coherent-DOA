clear; clc; close all;
% ============================================================
% FORMAL reduced ULA-transfer SNR comparison
% ST-BA vs CML vs FBSS-MUSIC
% ============================================================
thisFile=mfilename('fullpath'); root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project(); assert_project_dependencies(root); require_reduced_qualification(root);

cfg=config_v0(); cfg.K=2; cfg.p=0:9; cfg.M=10;
N=100; rho=0.8; snrList=-15:5:15; NMC=300;
fbOpts.K=2; fbOpts.thetaGrid=-60:0.1:60; fbOpts.subarrayLength=cfg.M-cfg.K+1;

methodNames={'ST-BA','CML','FBSS-MUSIC'}; nM=numel(methodNames); nS=numel(snrList);
mse=nan(NMC,nS,nM); mx=nan(NMC,nS,nM); runt=nan(NMC,nS,nM);

fprintf('============================================================\n');
fprintf(' FORMAL reduced ULA-transfer SNR comparison\n');
fprintf(' NMC=%d, N=%d, rho=%.2f, p=0:9\n',NMC,N,rho);
fprintf('============================================================\n');
for is=1:nS
    fprintf('SNR=%+.1f dB\n',snrList(is));
    for mc=1:NMC
        % Reuse the pre-specified V18 seed family; no SCCM formal outcomes
        % were generated before the method-set reduction.
        seed=cfg.seed+8400000+mc;
        [X,tr]=generate_trial_ar1_k(cfg,N,snrList(is),rho,[],seed);
        th=cell(1,nM);
        t0=tic; th{1}=estimate_stba_k2_runtime(X,cfg); runt(mc,is,1)=toc(t0);
        t0=tic; th{2}=cml_equalamp_k2_estimate(X,cfg.p); runt(mc,is,2)=toc(t0);
        t0=tic; th{3}=fbss_music_ula_estimate(X,cfg.p,fbOpts); runt(mc,is,3)=toc(t0);
        for im=1:nM
            [r,m]=angle_error_k(th{im},tr.thetaDeg); mse(mc,is,im)=r^2; mx(mc,is,im)=m;
        end
        if mod(mc,50)==0, fprintf('  completed %d/%d\n',mc,NMC); end
    end
end

[rmse,p05,med,p95,p99]=summarize(mse,mx);
ciSTminus=nan(nS,nM-1,2); Bboot=5000;
fprintf('\n==================== V8F RESULTS ====================\n');
for is=1:nS
    fprintf('\nSNR=%+.1f dB\n',snrList(is));
    for im=1:nM
        fprintf('  %-11s RMSE=%8.4f P0.5=%6.3f P95=%8.4f P99=%8.4f\n', ...
            methodNames{im},rmse(is,im),p05(is,im),p95(is,im),p99(is,im));
    end
    for im=2:nM
        d=mse(:,is,1)-mse(:,is,im);
        [lo,hi]=paired_bootstrap_ci(d,Bboot,21840000+100*is+im);
        ciSTminus(is,im-1,:)=[lo hi];
    end
end

results.version='V8F'; results.scope='candidate formal reduced ULA geometry-transfer SNR comparison';
results.cfg=cfg; results.fbssOptions=fbOpts; results.NMC=NMC; results.N=N; results.rho=rho; results.snrList=snrList;
results.methodNames=methodNames; results.doaMse=mse; results.maxErr=mx; results.runtime=runt;
results.rmse=rmse; results.p05=p05; results.medianMax=med; results.p95=p95; results.p99=p99;
results.ciSTminusBaseline_MSE=ciSTminus;
results.seedRule='cfg.seed + 8400000 + mc, paired across SNR values and methods';
results.protocolNote=['All methods receive the same 10-element ULA data, equal-magnitude fully coherent K=2 model, ' ...
    'and known K=2. CML uses the source-paper coarse/fine likelihood search. FBSS uses L=M-K+1.'];
results.excludedMethods={'SCCM-Schur-adapt','Toeplitz diagnostic'};
outDir=fullfile(root,'results','candidate'); if ~exist(outDir,'dir'),mkdir(outDir);end
save(fullfile(outDir,'v8f_ula_reduced_snr_results.mat'),'results','-v7.3');
fprintf('\nSaved candidate result. Existing frozen files were not modified.\n');

function require_reduced_qualification(root)
Q=load(fullfile(root,'results','qualification','v8e_ula_reduced_sanity_results.mat'),'results');
assert(Q.results.pass,'Reduced ULA sanity must pass first.');
end
function [rmse,p05,med,p95,p99]=summarize(mse,mx)
rmse=squeeze(sqrt(mean(mse,1))); p05=squeeze(mean(mx<=0.5,1)); med=squeeze(median(mx,1));
p95=squeeze(prctile(mx,95,1)); p99=squeeze(prctile(mx,99,1));
end
