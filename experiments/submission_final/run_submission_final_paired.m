function R = run_submission_final_paired(mode,projectRoot)
%RUN_SUBMISSION_FINAL_PAIRED Final fresh same-X evidence generation.
%
% Usage:
%   R=run_submission_final_paired('snr',pwd);
%   R=run_submission_final_paired('snapshot',pwd);
%   R=run_submission_final_paired('separation',pwd);
%   R=run_submission_final_paired('keypoint',pwd);
%   R=run_submission_final_paired('runtime',pwd);
%
% Recommended order:
%   snr -> snapshot -> separation -> keypoint -> runtime
%
% MAIN COMPARATORS
%   1) frozen original ST-BA
%   2) qualified continuous locally-polished equal-amplitude CML
%   3) equation-based TMR-VAI specialization
%
% EXCLUDED FROM MAIN CURVES
%   Peng independent reproduction: frozen not qualified
%   MESA independent reproduction: Phase-M2 frozen not qualified
%
% Core sweeps use NMC=1000.  The key point uses NMC=5000 specifically to
% audit heavy-tail Monte-Carlo stability.  All methods see exactly the same
% X within each trial.  No method receives truth DOAs.

if nargin<1 || isempty(mode),mode='snr';end
if nargin<2 || isempty(projectRoot),projectRoot=pwd;end
mode=lower(char(mode));projectRoot=char(projectRoot);

% Clean-package path setup.  The qualified continuous-CML wrapper lives
% in src/cml; no historical Peng-development folder is required.
addpath(projectRoot,'-begin');
setup_project();
thisDir=fileparts(mfilename('fullpath'));
addpath(thisDir,'-begin');
verify_submission_final_patch(projectRoot);

outDir=fullfile(projectRoot,'results','submission_final');
if exist(outDir,'dir')~=7,mkdir(outDir);end

cfg=config_v0();cfg.K=2;
methodNames={'ST-BA','Continuous CML','TMR-VAI'};
nMeth=3;

switch mode
    case 'snr'
        NMC=1000; N=100; rho=0.8;
        axisVals=[-15 -10 -5 0 5 10 15];
        nCond=numel(axisVals);
        seedBase=202621000;
        label='snr';
        gen=@(c,mc) generate_trial_ar1_k( ...
            cfg,N,axisVals(c),rho,[],seedBase+mc);

    case 'snapshot'
        NMC=1000; snrDb=-10; rho=0.8;
        axisVals=[50 100 200 300 500];
        nCond=numel(axisVals);
        seedBase=202622000;
        label='snapshot';
        % Same seed across N preserves angle/coherent-phase pairing.  Each
        % condition is still evaluated only by within-condition same-X
        % method pairing.
        gen=@(c,mc) generate_trial_ar1_k( ...
            cfg,axisVals(c),snrDb,rho,[],seedBase+mc);

    case 'separation'
        NMC=1000; N=100; snrDb=-10; rho=0.8;
        axisVals=[1 2 3 4 5 7.5 10];
        nCond=numel(axisVals);
        centerSeedBase=202623000;
        dataSeedBase=202624000;
        label='separation';
        gen=@(c,mc) local_sep_trial(cfg,N,snrDb,rho, ...
            axisVals(c),centerSeedBase+mc,dataSeedBase+mc);

    case 'keypoint'
        NMC=5000; N=100; snrDb=-10; rho=0.8;
        axisVals=1;nCond=1;seedBase=202625000;
        label='keypoint';
        gen=@(c,mc) generate_trial_ar1_k( ...
            cfg,N,snrDb,rho,[],seedBase+mc); %#ok<NASGU>

    case 'runtime'
        R=local_runtime(projectRoot,cfg,methodNames);
        return

    otherwise
        error('Unknown mode: %s',mode);
end

base=fullfile(outDir,['submission_final_' label]);
assert(exist([base '.mat'],'file')~=2, ...
    'Final result already exists; do not overwrite: %s.mat',base);

doaMse=nan(NMC,nCond,nMeth);
maxErr=nan(NMC,nCond,nMeth);
failed=false(NMC,nCond,nMeth);
runtimeSec=nan(NMC,nCond,nMeth);
thetaTrue=nan(NMC,nCond,2);
thetaHat=nan(NMC,nCond,nMeth,2);
crbMse=nan(NMC,nCond);
cmlAccepted=false(NMC,nCond);
cmlScoreGain=nan(NMC,nCond);

fprintf('\nSUBMISSION-FINAL %s PAIRED MONTE CARLO\n',upper(label));
fprintf('Methods: ST-BA / continuous CML / TMR-VAI\n');
fprintf('NMC=%d, conditions=%s\n',NMC,mat2str(axisVals));
fprintf('MESA and Peng remain excluded from main curves by frozen qualification status.\n');

for c=1:nCond
    fprintf('\ncondition %d/%d: %g\n',c,nCond,axisVals(c));
    for mc=1:NMC
        [X,truth]=gen(c,mc);
        rec=submission_final_eval_threeway(X,truth,cfg);

        thetaTrue(mc,c,:)=truth.thetaDeg;
        for m=1:nMeth
            doaMse(mc,c,m)=rec.doaMse(m);
            maxErr(mc,c,m)=rec.maxErr(m);
            failed(mc,c,m)=rec.failed(m);
            runtimeSec(mc,c,m)=rec.runtimeSec(m);
            thetaHat(mc,c,m,:)=rec.thetaHat(m,:);
        end
        crbMse(mc,c)=rec.crbMse;
        cmlAccepted(mc,c)=rec.cmlAccepted;
        cmlScoreGain(mc,c)=rec.cmlScoreGain;

        if mod(mc,50)==0
            fprintf('  completed %d/%d\n',mc,NMC);
        end
    end

    checkpoint=struct();
    checkpoint.label=label;checkpoint.axisVals=axisVals;
    checkpoint.completedCondition=c;
    checkpoint.NMC=NMC;checkpoint.methodNames=methodNames;
    checkpoint.doaMse=doaMse;checkpoint.maxErr=maxErr;
    checkpoint.failed=failed;checkpoint.runtimeSec=runtimeSec;
    checkpoint.thetaTrue=thetaTrue;checkpoint.thetaHat=thetaHat;
    checkpoint.crbMse=crbMse;checkpoint.cmlAccepted=cmlAccepted;
    checkpoint.cmlScoreGain=cmlScoreGain;
    save([base '_checkpoint.mat'],'checkpoint','-v7.3');
end

summary=submission_final_summarize(doaMse,maxErr,failed);
crbRmse=sqrt(mean(crbMse,1,'omitnan'));

R=struct();
R.version='submission-final-paired-v1';
R.label=label;
R.freshEvidence=true;
R.sameXWithinTrial=true;
R.NMC=NMC;
R.axisVals=axisVals;
R.methodNames=methodNames;
R.doaMse=doaMse;
R.maxErr=maxErr;
R.failed=failed;
R.runtimeSec=runtimeSec;
R.thetaTrue=thetaTrue;
R.thetaHat=thetaHat;
R.crbMse=crbMse;
R.crbRmse=crbRmse;
R.cmlAccepted=cmlAccepted;
R.cmlScoreGain=cmlScoreGain;
R.summary=summary;
R.protocol=local_protocol(mode);
R.excludedBaselines={ ...
    'Peng independent specialization: frozen not qualified', ...
    'MESA independent specialization: Phase-M2 frozen not qualified'};
R.posthocRetuning=false;

if strcmp(mode,'keypoint')
    R.mcConvergence=local_mc_convergence(doaMse,maxErr,failed);
end

save([base '.mat'],'R','-v7.3');
local_write_summary_csv(R,[base '_summary.csv']);

fprintf('\nSUMMARY (%s)\n',upper(label));
for c=1:nCond
    fprintf('condition=%g | CRB RMSE=%.5g\n',axisVals(c),crbRmse(c));
    for m=1:nMeth
        fprintf(['  %-15s RMSE=%8.4f P0.5=%6.3f P0.25=%6.3f ', ...
            'P95=%8.4f P99=%8.4f cat>5=%6.3f fail=%6.3f\n'], ...
            methodNames{m},summary.rmse(c,m),summary.p05(c,m), ...
            summary.p025(c,m),summary.p95Max(c,m),summary.p99Max(c,m), ...
            summary.catGt5(c,m),summary.failureRate(c,m));
    end
end
fprintf('Saved: %s.mat\n',base);
fprintf('Saved: %s_summary.csv\n',base);
end

function [X,truth]=local_sep_trial(cfg,N,snrDb,rho,sep,centerSeed,dataSeed)
old=rng;
cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(centerSeed,'twister');
center=-53+106*rand;
theta=[center-sep/2 center+sep/2];
[X,truth]=generate_trial_ar1_k(cfg,N,snrDb,rho,theta,dataSeed);
end

function P=local_protocol(mode)
P=struct();
P.mode=mode;
P.primaryArray=[0 1 2 3 4 5 11 17 23 29];
P.K=2;
P.mainRho=0.8;
P.equalMagnitudeCoherent=true;
P.mainMethods={'ST-BA','Continuous CML','TMR-VAI'};
P.coreNMC=1000;
P.keypointNMC=5000;
P.thresholdsDeg=[0.25 0.5 5];
P.tailPercentiles=[90 95 99];
P.noTruthInputToEstimators=true;
end

function C=local_mc_convergence(doaMse,maxErr,failed)
prefix=[100 300 500 1000 2000 5000];
nMeth=size(doaMse,3);
C=struct();
C.prefix=prefix;
C.rmse=nan(numel(prefix),nMeth);
C.p05=nan(numel(prefix),nMeth);
C.p99=nan(numel(prefix),nMeth);
C.catGt5=nan(numel(prefix),nMeth);
for i=1:numel(prefix)
    n=prefix(i);
    for m=1:nMeth
        ok=~failed(1:n,1,m) & isfinite(doaMse(1:n,1,m)) & isfinite(maxErr(1:n,1,m));
        if any(ok)
            C.rmse(i,m)=sqrt(mean(doaMse(ok,1,m)));
            x=maxErr(ok,1,m);
            C.p05(i,m)=mean(x<=0.5);
            C.p99(i,m)=prctile(x,99);
            C.catGt5(i,m)=mean(x>5);
        end
    end
end
end

function local_write_summary_csv(R,path)
S=R.summary;
rows=[];
for c=1:numel(R.axisVals)
    for m=1:numel(R.methodNames)
        rows=[rows; table(R.axisVals(c),string(R.methodNames{m}), ...
            S.rmse(c,m),S.p025(c,m),S.p05(c,m),S.medianMax(c,m), ...
            S.p90Max(c,m),S.p95Max(c,m),S.p99Max(c,m), ...
            S.catGt5(c,m),S.failureRate(c,m), ...
            'VariableNames',{'condition','method','rmseDeg','P025','P05', ...
            'medianMaxErrDeg','P90MaxErrDeg','P95MaxErrDeg','P99MaxErrDeg', ...
            'catastrophicGt5Rate','failureRate'})]; %#ok<AGROW>
    end
end
writetable(rows,path);
end

function R=local_runtime(projectRoot,cfg,methodNames)
outDir=fullfile(projectRoot,'results','submission_final');
if exist(outDir,'dir')~=7,mkdir(outDir);end
base=fullfile(outDir,'submission_final_runtime');
assert(exist([base '.mat'],'file')~=2,'Runtime output already exists.');

Ntrial=100;N=100;snrDb=-10;rho=0.8;seedBase=202626000;
times=nan(Ntrial,3);
failed=false(Ntrial,3);

% Warm-up: three fresh datasets, all methods once.
for j=1:3
    [X,tr]=generate_trial_ar1_k(cfg,N,snrDb,rho,[],seedBase-100+j);
    submission_final_eval_threeway(X,tr,cfg);
end

permList=[1 2 3;2 3 1;3 1 2;1 3 2;3 2 1;2 1 3];
sourceGrid=cfg.angleMin:0.5:cfg.angleMax;
musicGrid=cfg.angleMin:0.1:cfg.angleMax;

fprintf('\nSUBMISSION-FINAL MATCHED RUNTIME\n');
for i=1:Ntrial
    [X,~]=generate_trial_ar1_k(cfg,N,snrDb,rho,[],seedBase+i);
    ord=permList(mod(i-1,size(permList,1))+1,:);
    for jj=1:3
        m=ord(jj);
        try
            t=tic;
            switch m
                case 1
                    estimate_stba_k2_runtime(X,cfg);
                case 2
                    cml_equalamp_k2_continuous_polish(X,cfg.p);
                case 3
                    tmr_vai_fully_augmented_estimate( ...
                        X,cfg.p,2,sourceGrid,musicGrid,'mean');
            end
            times(i,m)=toc(t);
        catch
            failed(i,m)=true;
        end
    end
    if mod(i,20)==0,fprintf('  completed %d/%d\n',i,Ntrial);end
end

medianSec=median(times,1,'omitnan');
p95Sec=prctile(times,95,1);

R=struct();
R.version='submission-final-runtime-v1';
R.Ntrial=Ntrial;R.N=N;R.snrDb=snrDb;R.rho=rho;
R.methodNames=methodNames;R.timesSec=times;R.failed=failed;
R.medianSec=medianSec;R.p95Sec=p95Sec;
R.stOverCml=medianSec(1)/medianSec(2);
R.stOverTmr=medianSec(1)/medianSec(3);
R.cmlOverTmr=medianSec(2)/medianSec(3);

save([base '.mat'],'R','-v7.3');
T=table(string(methodNames(:)),medianSec(:),p95Sec(:),mean(failed,1).', ...
    'VariableNames',{'method','medianSec','p95Sec','failureRate'});
writetable(T,[base '.csv']);

fprintf('  ST-BA median %.6fs\n',medianSec(1));
fprintf('  Continuous CML median %.6fs\n',medianSec(2));
fprintf('  TMR-VAI median %.6fs\n',medianSec(3));
fprintf('  ST/CML=%.4f ST/TMR=%.4f CML/TMR=%.4f\n', ...
    R.stOverCml,R.stOverTmr,R.cmlOverTmr);
fprintf('Saved: %s.mat/.csv\n',base);
end
