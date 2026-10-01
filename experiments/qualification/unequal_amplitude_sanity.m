clear; clc; close all;
% Unequal-amplitude ST-BA qualification. No performance-ranking gate.
thisFile=mfilename('fullpath'); root=fileparts(fileparts(fileparts(thisFile)));
addpath(root,'-begin'); setup_project(); assert_project_dependencies(root);
cfg=config_v0(); cfg.K=2;
N=100; snrDb=-10; rho=0.8; ratioDb=[0 -6 -9]; NMC=50;

% Generator compatibility: equal magnitudes must reproduce the established generator.
seed0=cfg.seed+9099001;
[Xa,Ta]=generate_trial_ar1_k(cfg,N,snrDb,rho,[],seed0);
[Xb,Tb]=generate_trial_ar1_k_weighted(cfg,N,snrDb,rho,[],seed0,[1;1]);
compatX=norm(Xa-Xb,'fro') <= 1e-12*max(1,norm(Xa,'fro'));
compatTheta=isequaln(Ta.thetaDeg,Tb.thetaDeg);
compatBeta=norm(Ta.beta-Tb.beta) <= 1e-13;
assert(compatX && compatTheta && compatBeta,'Weighted generator does not reproduce equal-magnitude control.');

mse=nan(NMC,numel(ratioDb)); mx=mse; bnmse=mse; sel=nan(NMC,numel(ratioDb)); kMargin=sel;
finiteOK=true;
fprintf('============================================================\n');
fprintf(' Unequal-amplitude sanity | NMC=%d N=%d SNR=%g dB rho=%.2f\n',NMC,N,snrDb,rho);
fprintf(' ratio dB = [%s]\n',num2str(ratioDb));
fprintf('============================================================\n');
for ir=1:numel(ratioDb)
    r=10^(ratioDb(ir)/20);
    for mc=1:NMC
        % Balance which ordered source is weaker; no extra RNG draw.
        weakIndex=1+mod(mc,2); betaMag=ones(2,1); betaMag(weakIndex)=r;
        seed=cfg.seed+9000000+mc; % paired across amplitude-ratio conditions
        [X,tr]=generate_trial_ar1_k_weighted(cfg,N,snrDb,rho,[],seed,betaMag);
        try
            [th,info]=estimate_stba_k2_runtime(X,cfg);
            [er,em]=angle_error_k(th,tr.thetaDeg);
            mse(mc,ir)=er^2; mx(mc,ir)=em; bnmse(mc,ir)=10*log10(max(nmse_B(info.B,tr.Btrue),realmin));
            sel(mc,ir)=info.stability.selectedIndex; kMargin(mc,ir)=info.stability.kappaMargin;
        catch ME
            finiteOK=false; warning('Unequal-amplitude sanity failure ratio=%g dB mc=%d: %s',ratioDb(ir),mc,ME.message);
        end
    end
end
rmse=sqrt(mean(mse,1)); p05=mean(mx<=0.5,1); zlpRate=mean(sel==2,1);
pass=finiteOK && all(isfinite(mse(:))) && all(isfinite(mx(:))) && all(isfinite(bnmse(:))) && all(ismember(sel(:),[1 2]));
fprintf('\nratio(dB)    RMSE(deg)   P0.5   B-NMSE(dB)  ZL-select\n');
for ir=1:numel(ratioDb)
    fprintf('%8.1f    %9.4f  %6.3f   %10.3f   %8.3f\n',ratioDb(ir),rmse(ir),p05(ir),mean(bnmse(:,ir)),zlpRate(ir));
end
results.version='UA-Q1'; results.pass=pass; results.scope='qualification-only unequal-amplitude ST-BA sweep';
results.cfg=cfg; results.N=N; results.snrDb=snrDb; results.rho=rho; results.ratioDb=ratioDb; results.NMC=NMC;
results.doaMse=mse; results.maxErr=mx; results.bNmseDb=bnmse; results.selectedIndex=sel; results.kappaMargin=kMargin;
results.rmse=rmse; results.p05=p05; results.zlpRate=zlpRate;
results.generatorCompatibility=struct('X',compatX,'theta',compatTheta,'beta',compatBeta);
results.seedRule='cfg.seed + 9000000 + mc; same latent trial reused across amplitude-ratio conditions; weak source alternates by mc parity';
results.performanceGate='NONE. PASS checks generator compatibility and finite ST-BA outputs only.';
outDir=fullfile(root,'results','qualification'); save(fullfile(outDir,'unequal_amplitude_sanity_results.mat'),'results','-v7.3');
if pass, fprintf('\nPASS: unequal-amplitude pipeline is eligible for formal characterization.\n'); else, error('Unequal-amplitude sanity failed.'); end
