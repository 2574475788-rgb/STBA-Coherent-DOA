function msg=verify_submission_final_patch(projectRoot)
%VERIFY_SUBMISSION_FINAL_PATCH
% Preflight for the final fresh paired Monte-Carlo evidence generation.
%
% This stage does not modify ST-BA, CML, TMR-VAI, MESA or Peng algorithms.

if nargin<1 || isempty(projectRoot),projectRoot=pwd;end
projectRoot=char(projectRoot);

reqFns={ ...
    'generate_trial_ar1_k', ...
    'estimate_stba_k2_runtime', ...
    'cml_equalamp_k2_continuous_polish', ...
    'tmr_vai_fully_augmented_estimate', ...
    'coherent_equalamp_constrained_crb_k2', ...
    'angle_error_k', ...
    'paired_bootstrap_ci', ...
    'submission_final_eval_threeway', ...
    'submission_final_summarize'};
for k=1:numel(reqFns)
    assert(exist(reqFns{k},'file')==2, ...
        'Required function missing from path: %s',reqFns{k});
end

cfg=config_v0();
assert(cfg.K==2);
assert(isequal(cfg.p,[0 1 2 3 4 5 11 17 23 29]));

% Ensure the clean repository resolves the qualified continuous-polish
% wrapper from src/cml.
cmlPath=which('cml_equalamp_k2_continuous_polish');
expectedCml=fullfile(projectRoot,'src','cml','cml_equalamp_k2_continuous_polish.m');
assert(strcmpi(strrep(cmlPath,'/','\'),strrep(expectedCml,'/','\')), ...
    'Unexpected continuous CML implementation on MATLAB path: %s',cmlPath);

% Retained fixed supporting evidence used by publication figures.
assert(exist(fullfile(projectRoot,'results','frozen','v1g_stability_selector_results.mat'),'file')==2);
assert(exist(fullfile(projectRoot,'results','frozen','v3c_robustness_formal_results.mat'),'file')==2);

msg=['SUBMISSION-FINAL PREFLIGHT PASS; algorithms are frozen. ', ...
    'Fresh paired evidence will be written only to results/submission_final.'];
fprintf('%s\n',msg);
end
