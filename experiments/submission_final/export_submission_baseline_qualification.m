function T = export_submission_baseline_qualification(projectRoot)
%EXPORT_SUBMISSION_BASELINE_QUALIFICATION
% Create an auditable supplementary table for external baseline qualification.
% This does not turn failed independent reproductions into paper baselines.

if nargin<1 || isempty(projectRoot),projectRoot=pwd;end
projectRoot=char(projectRoot);
outDir=fullfile(projectRoot,'results','submission_final');
if exist(outDir,'dir')~=7,mkdir(outDir);end

name=["Peng PCA-VP independent"; ...
      "Peng MANM independent development"; ...
      "MESA independent"];
status=["NOT QUALIFIED";"NOT QUALIFIED";"NOT QUALIFIED"];
authorCode=[false;false;false];
mainCurve=[false;false;false];
evidence=["Phase-2W frozen 23/24 numerical closure; not published-Peng baseline"; ...
          "Phase-3D/3E development failed with catastrophic cases"; ...
          "Phase-M2: 12/12 numerical protocol, 4/4 paper-inspired <=2deg, project 6/8 <=2deg, one >5deg"];
notes=["May be discussed only as independent PCA-VP specialization"; ...
       "Post-hoc dual-polynomial audit did not create qualification"; ...
       "Do not retune mu from Phase-M2 outcomes"];

T=table(name,status,authorCode,mainCurve,evidence,notes);
writetable(T,fullfile(outDir,'external_baseline_qualification.csv'));

fid=fopen(fullfile(outDir,'external_baseline_qualification_README.txt'),'w');
fprintf(fid,['External baseline qualification record\n\n', ...
    'Peng and MESA were not silently omitted after performance comparison. ', ...
    'Their independent implementations were subjected to pre-specified ', ...
    'qualification protocols and frozen when those protocols failed. ', ...
    'They are therefore not promoted to validated main performance curves.\n']);
fclose(fid);
fprintf('Saved external baseline qualification table in results/submission_final.\n');
end
