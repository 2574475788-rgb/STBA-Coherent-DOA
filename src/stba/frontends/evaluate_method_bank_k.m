function metrics = evaluate_method_bank_k(bank, truth, cfg, K, L, beamWidth, minStartDist)
%EVALUATE_METHOD_BANK_K Evaluate all frozen B estimates with the SAME VP backend.

if nargin < 5 || isempty(L), L=8; end
if nargin < 6 || isempty(beamWidth), beamWidth=80; end
if nargin < 7 || isempty(minStartDist), minStartDist=1.5; end

nMethods = numel(bank.names);

metrics.nmse = nan(1,nMethods);
metrics.doaMse = nan(1,nMethods);
metrics.maxErr = nan(1,nMethods);
metrics.vpResidual = nan(1,nMethods);
metrics.thetaHat = nan(nMethods,K);

for m=1:nMethods
    B = bank.B{m};

    metrics.nmse(m) = nmse_B(B,truth.Btrue);

    [th,info] = estimate_doa_from_B_varpro_k( ...
        B,cfg,K,L,beamWidth,minStartDist);

    [r,mx] = angle_error_k(th,truth.thetaDeg);

    metrics.doaMse(m) = r^2;
    metrics.maxErr(m) = mx;
    metrics.vpResidual(m) = info.finalResidual;
    metrics.thetaHat(m,:) = th;
end
end
