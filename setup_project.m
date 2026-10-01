function root = setup_project()
%SETUP_PROJECT Add only the curated reproducibility folders to MATLAB path.
%
% This clean repository intentionally excludes historical patch, diagnostic,
% failed-baseline-development, checkpoint and third-party folders.

root=fileparts(mfilename('fullpath'));

addpath(root,'-begin');
addpath(genpath(fullfile(root,'src')),'-begin');
addpath(fullfile(root,'experiments','submission_final'),'-begin');
addpath(fullfile(root,'experiments','formal'),'-begin');
addpath(fullfile(root,'experiments','qualification'),'-begin');
addpath(fullfile(root,'experiments','reproduction'),'-begin');
rehash path;

fprintf('ML_CFR clean reproducibility paths added from:\n  %s\n',root);
end
