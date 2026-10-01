clear; clc; close all;
% Reduced ULA geometry-transfer sanity: ST-BA / CML / FBSS-MUSIC.
% Qualification only. PASS is finite-output validity, not performance ranking.
thisFile=mfilename('fullpath'); root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project(); assert_project_dependencies(root);
cfg=config_v0(); cfg.K=2; cfg.p=0:9; cfg.M=10;
N=100; rho=0.8; snrList=[-10 0 10]; NMC=50;
fbOpts.K=2; fbOpts.thetaGrid=-60:0.1:60; fbOpts.subarrayLength=cfg.M-cfg.K+1;
methodNames={'ST-BA','CML','FBSS-MUSIC'}; nM=numel(methodNames); nS=numel(snrList);
mse=nan(NMC,nS,nM); mx=nan(NMC,nS,nM); runt=nan(NMC,nS,nM); finiteOK=true;
fprintf('============================================================\n');
fprintf(' Reduced ULA sanity: %d trials/condition\n',NMC);
fprintf(' Methods: ST-BA / CML / FBSS-MUSIC\n');
fprintf(' N=%d, rho=%.2f, p=0:9, FBSS L=%d Q=%d\n',N,rho,fbOpts.subarrayLength,cfg.M-fbOpts.subarrayLength+1);
fprintf('============================================================\n');
for is=1:nS
    fprintf('SNR=%+.1f dB\n',snrList(is));
    for mc=1:NMC
        seed=cfg.seed+8800000+100000*is+mc;
        [X,tr]=generate_trial_ar1_k(cfg,N,snrList(is),rho,[],seed);
        th=cell(1,nM);
        try, t0=tic; th{1}=estimate_stba_k2_runtime(X,cfg); runt(mc,is,1)=toc(t0); catch ME, finiteOK=false; warning('ST-BA failure SNR=%g mc=%d: %s',snrList(is),mc,ME.message); end
        try, t0=tic; th{2}=cml_equalamp_k2_estimate(X,cfg.p); runt(mc,is,2)=toc(t0); catch ME, finiteOK=false; warning('CML failure SNR=%g mc=%d: %s',snrList(is),mc,ME.message); end
        try, t0=tic; th{3}=fbss_music_ula_estimate(X,cfg.p,fbOpts); runt(mc,is,3)=toc(t0); catch ME, finiteOK=false; warning('FBSS failure SNR=%g mc=%d: %s',snrList(is),mc,ME.message); end
        for im=1:nM
            if isempty(th{im}) || numel(th{im})~=2 || any(~isfinite(th{im})), finiteOK=false; continue; end
            [r,m]=angle_error_k(th{im},tr.thetaDeg); mse(mc,is,im)=r^2; mx(mc,is,im)=m;
        end
    end
end
rmse=squeeze(sqrt(mean(mse,1,'omitnan'))); p05=squeeze(mean(mx<=0.5,1,'omitnan'));
med=squeeze(median(mx,1,'omitnan')); p95=squeeze(prctile(mx,95,1)); p99=squeeze(prctile(mx,99,1));
fprintf('\n==================== RESULTS ====================\n');
for is=1:nS
    fprintf('\nSNR=%+.1f dB\n',snrList(is));
    for im=1:nM
        fprintf('  %-11s RMSE=%8.4f P0.5=%6.3f medMax=%8.4f P95=%8.4f P99=%8.4f time=%8.4fs\n',methodNames{im},rmse(is,im),p05(is,im),med(is,im),p95(is,im),p99(is,im),median(runt(:,is,im),'omitnan'));
    end
end
pass=finiteOK && all(isfinite(mse(:))) && all(isfinite(mx(:)));
results.version='ULA-Q-clean'; results.scope='qualification-only reduced ULA transfer sanity'; results.pass=pass;
results.cfg=cfg; results.N=N; results.rho=rho; results.snrList=snrList; results.NMC=NMC; results.methodNames=methodNames; results.fbssOptions=fbOpts;
results.doaMse=mse; results.maxErr=mx; results.runtime=runt; results.rmse=rmse; results.p05=p05; results.medianMax=med; results.p95=p95; results.p99=p99;
results.seedRule='cfg.seed + 8800000 + 100000*SNRindex + mc; methods paired within trial';
results.performanceGate='NONE. PASS requires finite, valid outputs on every sanity trial.';
results.methodDecision='SCCM/Toeplitz diagnostic code removed after fidelity failure; active ULA set is ST-BA/CML/FBSS-MUSIC.';
outDir=fullfile(root,'results','qualification'); if ~exist(outDir,'dir'),mkdir(outDir);end
save(fullfile(outDir,'v8e_ula_reduced_sanity_results.mat'),'results','-v7.3');
if pass, fprintf('\nPASS: all three active methods produced finite outputs.\n'); else, error('Reduced ULA sanity failed.'); end
