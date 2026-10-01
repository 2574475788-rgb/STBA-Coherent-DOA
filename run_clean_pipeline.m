function run_clean_pipeline(action)
%RUN_CLEAN_PIPELINE Entry point for the public GitHub/Zenodo package.
%
% Actions
%   'check'       dependency + required-result integrity checks
%   'postprocess' recompute paired bootstrap tables from included final
%                 Monte-Carlo records; no estimator rerun
%   'core'        run the final paired Monte-Carlo modes only when the
%                 corresponding output file is absent
%
% This public package intentionally omits publication-figure layout scripts.
% The included MAT/CSV files are the numerical data underlying the paper.

if nargin<1 || isempty(action), action='check'; end
action=lower(char(action));

root=fileparts(mfilename('fullpath'));
cd(root);
setup_project();
assert_project_dependencies(root);

switch action
    case 'check'
        local_check_results(root);
        verify_submission_final_patch(root);
        fprintf('\nPUBLIC PACKAGE CHECK PASS.\n');

    case 'postprocess'
        analyze_submission_final_paired_ci(root);
        export_submission_baseline_qualification(root);
        fprintf('\nPostprocessing complete.\n');

    case 'core'
        modes={'snr','snapshot','separation','keypoint','runtime'};
        for k=1:numel(modes)
            out=fullfile(root,'results','submission_final', ...
                ['submission_final_' modes{k} '.mat']);
            if exist(out,'file')==2
                fprintf('SKIP %-10s existing result: %s\n',modes{k},out);
            else
                run_submission_final_paired(modes{k},root);
            end
        end

    otherwise
        error('Unknown action: %s. Use check, postprocess, or core.',action);
end
end

function local_check_results(root)
req={ ...
    fullfile(root,'results','submission_final','submission_final_snr.mat'); ...
    fullfile(root,'results','submission_final','submission_final_snapshot.mat'); ...
    fullfile(root,'results','submission_final','submission_final_separation.mat'); ...
    fullfile(root,'results','submission_final','submission_final_keypoint.mat'); ...
    fullfile(root,'results','submission_final','submission_final_runtime.mat'); ...
    fullfile(root,'results','frozen','v1g_stability_selector_results.mat'); ...
    fullfile(root,'results','frozen','v3c_robustness_formal_results.mat'); ...
    fullfile(root,'results','candidate','v8f_ula_reduced_snr_results.mat'); ...
    fullfile(root,'results','candidate','v8g_ula_reduced_snapshot_results.mat'); ...
    fullfile(root,'results','candidate','unequal_amplitude_formal_results.mat'); ...
    fullfile(root,'results','qualification','v8e_ula_reduced_sanity_results.mat'); ...
    fullfile(root,'results','qualification','unequal_amplitude_sanity_results.mat')};
for k=1:numel(req)
    assert(exist(req{k},'file')==2,'Missing required stored result: %s',req{k});
end
assert_ula_candidate_integrity(root);
fprintf('PASS: required stored result files are present.\n');
end
