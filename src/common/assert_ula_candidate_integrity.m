function report = assert_ula_candidate_integrity(root)
%ASSERT_ULA_CANDIDATE_INTEGRITY Verify reduced ULA candidate result set.
if nargin<1 || isempty(root), root=fileparts(fileparts(fileparts(mfilename('fullpath')))); end
Q=load(fullfile(root,'results','qualification','v8e_ula_reduced_sanity_results.mat'),'results');
assert(Q.results.pass,'Reduced ULA sanity qualification did not pass.');
A=load(fullfile(root,'results','candidate','v8f_ula_reduced_snr_results.mat'),'results'); A=A.results;
B=load(fullfile(root,'results','candidate','v8g_ula_reduced_snapshot_results.mat'),'results'); B=B.results;
expected={'ST-BA','CML','FBSS-MUSIC'};
checks.methods=isequal(A.methodNames,expected) && isequal(B.methodNames,expected);
checks.finite=all(isfinite(A.doaMse(:))) && all(isfinite(A.maxErr(:))) && all(isfinite(B.doaMse(:))) && all(isfinite(B.maxErr(:)));
checks.axes=isequal(A.snrList,-15:5:15) && isequal(B.NList,[50 100 200 300 500]);
checks.cardinality=size(A.doaMse,1)==300 && size(B.doaMse,1)==300;
report.checks=checks; report.pass=all(structfun(@logical,checks));
assert(report.pass,'Reduced ULA candidate integrity check failed.');
end
