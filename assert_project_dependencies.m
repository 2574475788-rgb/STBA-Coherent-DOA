function report = assert_project_dependencies(root)
%ASSERT_PROJECT_DEPENDENCIES Fail fast on missing/shadowed clean dependencies.

if nargin<1 || isempty(root)
    root=fileparts(mfilename('fullpath'));
end

req = {
    'config_v0',                         fullfile(root,'src','common','config_v0.m');
    'angle_error_k',                     fullfile(root,'src','common','angle_error_k.m');
    'nmse_B',                            fullfile(root,'src','common','nmse_B.m');
    'steering_matrix',                   fullfile(root,'src','common','steering_matrix.m');
    'generate_trial_ar1_k',              fullfile(root,'src','generators','generate_trial_ar1_k.m');
    'generate_trial_ar1_k_weighted',     fullfile(root,'src','generators','generate_trial_ar1_k_weighted.m');
    'generate_trial_mismatch_k',         fullfile(root,'src','generators','generate_trial_mismatch_k.m');
    'estimate_stba_k2_runtime',          fullfile(root,'src','stba','frontends','estimate_stba_k2_runtime.m');
    'stability_adaptive_ba_fgls',        fullfile(root,'src','stba','frontends','stability_adaptive_ba_fgls.m');
    'estimate_doa_from_B_varpro_k',      fullfile(root,'src','stba','backend','estimate_doa_from_B_varpro_k.m');
    'cml_equalamp_k2_estimate',          fullfile(root,'src','cml','cml_equalamp_k2_estimate.m');
    'cml_equalamp_k2_continuous_polish', fullfile(root,'src','cml','cml_equalamp_k2_continuous_polish.m');
    'tmr_vai_fully_augmented_estimate',  fullfile(root,'src','tmr','tmr_vai_fully_augmented_estimate.m');
    'fbss_music_ula_estimate',           fullfile(root,'src','baselines','fbss_music_ula_estimate.m');
    'coherent_equalamp_constrained_crb_k2', fullfile(root,'src','crb','coherent_equalamp_constrained_crb_k2.m');
    };

fprintf('Clean-project dependency preflight:\n');
resolved=cell(size(req,1),1);
pass=false(size(req,1),1);
for k=1:size(req,1)
    name=req{k,1}; expected=req{k,2};
    assert(exist(expected,'file')==2,'Missing clean-project file:\n%s',expected);
    got=which(name); resolved{k}=got;
    pass(k)=same_file(got,expected);
    fprintf('  %-38s [%s]\n',name,passfail(pass(k)));
end
assert(all(pass),'A required clean-project dependency is missing or shadowed.');
report=struct('pass',true,'names',{req(:,1)},'resolved',{resolved},'root',root);
fprintf('PASS: clean dependency set resolves only to this repository.\n');
end

function tf=same_file(a,b)
if isempty(a), tf=false; return; end
a=strrep(char(a),'/','\');
b=strrep(char(b),'/','\');
tf=strcmpi(a,b);
end

function s=passfail(tf)
if tf, s='PASS'; else, s='FAIL'; end
end
