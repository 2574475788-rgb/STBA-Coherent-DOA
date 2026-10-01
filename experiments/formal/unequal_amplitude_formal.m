clear; clc; close all;
% Formal ST-BA operating-boundary sweep versus coherent-source amplitude ratio.
thisFile=mfilename('fullpath'); root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project(); assert_project_dependencies(root);
Q=load(fullfile(root,'results','qualification','unequal_amplitude_sanity_results.mat'),'results');
assert(Q.results.pass,'Run unequal-amplitude qualification first.');
cfg=config_v0(); cfg.K=2;
N=100; snrDb=-10; rho=0.8; ratioDb=[0 -3 -6 -9]; ratioLin=10.^(ratioDb/20); NMC=300; Bboot=5000;
nR=numel(ratioDb); mse=nan(NMC,nR); mx=mse; bnmse=mse; sel=nan(NMC,nR); kMargin=sel;
fprintf('============================================================\n');
fprintf(' FORMAL unequal-amplitude ST-BA characterization\n');
fprintf(' NMC=%d N=%d SNR=%g dB rho=%.2f ratios=[%s] dB\n',NMC,N,snrDb,rho,num2str(ratioDb));
fprintf('============================================================\n');
for ir=1:nR
    fprintf('ratio=%+.1f dB\n',ratioDb(ir)); r=ratioLin(ir);
    for mc=1:NMC
        weakIndex=1+mod(mc,2); betaMag=ones(2,1); betaMag(weakIndex)=r;
        seed=cfg.seed+9100000+mc; % paired across ratio conditions
        [X,tr]=generate_trial_ar1_k_weighted(cfg,N,snrDb,rho,[],seed,betaMag);
        [th,info]=estimate_stba_k2_runtime(X,cfg);
        [er,em]=angle_error_k(th,tr.thetaDeg);
        mse(mc,ir)=er^2; mx(mc,ir)=em; bnmse(mc,ir)=10*log10(max(nmse_B(info.B,tr.Btrue),realmin));
        sel(mc,ir)=info.stability.selectedIndex; kMargin(mc,ir)=info.stability.kappaMargin;
        if mod(mc,50)==0, fprintf('  completed %d/%d\n',mc,NMC); end
    end
end
rmse=sqrt(mean(mse,1)); p05=mean(mx<=0.5,1); med=median(mx,1); p95=prctile(mx,95,1); p99=prctile(mx,99,1);
meanBNmse=mean(bnmse,1); zlpRate=mean(sel==2,1);
% Paired changes relative to equal amplitude (0 dB), with pre-specified bootstrap.
ciMseChange=nan(nR,2); ciSuccessChange=nan(nR,2);
baseSuccess=double(mx(:,1)<=0.5);
for ir=1:nR
    [ciMseChange(ir,1),ciMseChange(ir,2)]=paired_bootstrap_ci(mse(:,ir)-mse(:,1),Bboot,29100000+ir);
    [ciSuccessChange(ir,1),ciSuccessChange(ir,2)]=paired_bootstrap_ci(double(mx(:,ir)<=0.5)-baseSuccess,Bboot,29110000+ir);
end
fprintf('\n==================== RESULTS ====================\n');
fprintf('ratio(dB)  RMSE(deg)   P0.5  medMax  P95   P99   B-NMSE(dB)  ZL-select\n');
for ir=1:nR
    fprintf('%8.1f  %9.4f  %6.3f  %6.3f %6.3f %6.3f   %10.3f   %8.3f\n',ratioDb(ir),rmse(ir),p05(ir),med(ir),p95(ir),p99(ir),meanBNmse(ir),zlpRate(ir));
end
results.version='UA-F1'; results.scope='candidate formal ST-BA unequal-amplitude operating-boundary characterization';
results.cfg=cfg; results.N=N; results.snrDb=snrDb; results.rho=rho; results.NMC=NMC; results.Bboot=Bboot;
results.ratioDb=ratioDb; results.ratioLin=ratioLin; results.doaMse=mse; results.maxErr=mx; results.bNmseDb=bnmse; results.selectedIndex=sel; results.kappaMargin=kMargin;
results.rmse=rmse; results.p05=p05; results.medianMax=med; results.p95=p95; results.p99=p99; results.meanBNmseDb=meanBNmse; results.zlpRate=zlpRate;
results.ciMseChangeVsEqual=ciMseChange; results.ciSuccessChangeVsEqual=ciSuccessChange;
results.seedRule='cfg.seed + 9100000 + mc; exactly paired latent DOA/phase/waveform/noise direction across amplitude ratios; weak source alternates by mc parity';
results.metricNote='Total finite-record array-input SNR is held fixed separately at each amplitude ratio; the weaker source therefore loses relative contribution without changing total SNR.';
results.scopeNote='ST-BA scope qualification only. Equal-amplitude CML is intentionally not used under unequal amplitudes.';
outDir=fullfile(root,'results','candidate'); save(fullfile(outDir,'unequal_amplitude_formal_results.mat'),'results','-v7.3');
fprintf('\nSaved candidate unequal-amplitude result. Frozen results were not modified.\n');
