function A = analyze_submission_final_paired_ci(projectRoot)
%ANALYZE_SUBMISSION_FINAL_PAIRED_CI
% Post-process the already completed submission-final Monte Carlo.
% NO estimator is rerun and NO random observation is regenerated.
%
% Computes paired percentile-bootstrap confidence intervals for:
%   ST-BA - Continuous CML
% on:
%   RMSE difference
%   P0.5 difference
%   catastrophic (>5 deg) rate difference
%
% It also exports the keypoint prefix-convergence table.
%
% Sign convention for differences:
%   dRMSE < 0  => ST-BA lower RMSE
%   dP05  > 0  => ST-BA higher P0.5
%   dCat  < 0  => ST-BA lower catastrophic rate

if nargin<1 || isempty(projectRoot),projectRoot=pwd;end
projectRoot=char(projectRoot);
resDir=fullfile(projectRoot,'results','submission_final');
assert(exist(resDir,'dir')==7,'Missing results/submission_final.');

files={ ...
    'submission_final_snr.mat','snr'; ...
    'submission_final_snapshot.mat','snapshot'; ...
    'submission_final_separation.mat','separation'; ...
    'submission_final_keypoint.mat','keypoint'};

B=10000;
rng(20260930,'twister');
rows=[];

for ff=1:size(files,1)
    D=load(fullfile(resDir,files{ff,1}),'R');
    R=D.R;
    for c=1:numel(R.axisVals)
        mse=squeeze(R.doaMse(:,c,:));
        mx=squeeze(R.maxErr(:,c,:));
        assert(size(mse,2)>=2 && size(mx,2)>=2);

        dRmse=sqrt(mean(mse(:,1),'omitnan'))-sqrt(mean(mse(:,2),'omitnan'));
        dP05=mean(mx(:,1)<=0.5)-mean(mx(:,2)<=0.5);
        dCat=mean(mx(:,1)>5)-mean(mx(:,2)>5);

        n=size(mse,1);
        br=nan(B,1);bp=nan(B,1);bc=nan(B,1);
        for b=1:B
            ii=randi(n,n,1);
            br(b)=sqrt(mean(mse(ii,1),'omitnan'))-sqrt(mean(mse(ii,2),'omitnan'));
            bp(b)=mean(mx(ii,1)<=0.5)-mean(mx(ii,2)<=0.5);
            bc(b)=mean(mx(ii,1)>5)-mean(mx(ii,2)>5);
        end

        r=table(string(files{ff,2}),R.axisVals(c),n, ...
            dRmse,prctile(br,2.5),prctile(br,97.5), ...
            dP05,prctile(bp,2.5),prctile(bp,97.5), ...
            dCat,prctile(bc,2.5),prctile(bc,97.5), ...
            'VariableNames',{'family','condition','NMC', ...
            'dRMSE_STminusCML','dRMSE_CIlo','dRMSE_CIhi', ...
            'dP05_STminusCML','dP05_CIlo','dP05_CIhi', ...
            'dCatGt5_STminusCML','dCat_CIlo','dCat_CIhi'});
        rows=[rows;r]; %#ok<AGROW>
    end
end

ciPath=fullfile(resDir,'submission_final_paired_STBA_vs_CML_CI.csv');
writetable(rows,ciPath);

K=load(fullfile(resDir,'submission_final_keypoint.mat'),'R');
C=K.R.mcConvergence;
conv=table(C.prefix(:),C.rmse(:,1),C.rmse(:,2),C.rmse(:,3), ...
    C.p05(:,1),C.p05(:,2),C.p05(:,3), ...
    C.p99(:,1),C.p99(:,2),C.p99(:,3), ...
    C.catGt5(:,1),C.catGt5(:,2),C.catGt5(:,3), ...
    'VariableNames',{'NMC', ...
    'RMSE_STBA','RMSE_CML','RMSE_TMR', ...
    'P05_STBA','P05_CML','P05_TMR', ...
    'P99_STBA','P99_CML','P99_TMR', ...
    'CatGt5_STBA','CatGt5_CML','CatGt5_TMR'});
convPath=fullfile(resDir,'submission_final_keypoint_MC_convergence.csv');
writetable(conv,convPath);

A=struct();
A.version='submission-final-postprocess-v1';
A.noEstimatorRerun=true;
A.noNewMonteCarlo=true;
A.bootstrapReplicates=B;
A.pairedCI=rows;
A.keypointConvergence=conv;
A.status='POSTPROCESS COMPLETE; READY FOR MANUSCRIPT/FIGURE REWRITE';
save(fullfile(resDir,'submission_final_postprocess.mat'),'A','-v7.3');

fprintf('\nSUBMISSION-FINAL POSTPROCESS\n');
disp(rows);
fprintf('\nKEYPOINT MONTE-CARLO CONVERGENCE\n');
disp(conv);
fprintf('Saved: %s\n',ciPath);
fprintf('Saved: %s\n',convPath);
fprintf('STATUS: %s\n',A.status);
end
