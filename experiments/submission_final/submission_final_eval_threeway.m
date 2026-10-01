function rec = submission_final_eval_threeway(X,truth,cfg)
%SUBMISSION_FINAL_EVAL_THREEWAY
% Same-X evaluation of the three qualified/retained paper comparators:
%   ST-BA (frozen original method)
%   continuous locally-polished equal-amplitude CML
%   equation-based TMR-VAI specialization
%
% MESA and Peng are intentionally NOT evaluated here:
% their independent implementations failed their own pre-fixed qualification
% protocols and are archived as qualification evidence, not promoted to
% validated performance curves.

sourceGrid=cfg.angleMin:0.5:cfg.angleMax;
musicGrid=cfg.angleMin:0.1:cfg.angleMax;

rec=struct();
rec.thetaTrue=truth.thetaDeg(:).';
rec.methods={'STBA','CML_CONTINUOUS','TMR_VAI'};
rec.thetaHat=nan(3,2);
rec.doaMse=nan(3,1);
rec.maxErr=nan(3,1);
rec.runtimeSec=nan(3,1);
rec.failed=false(3,1);
rec.status=strings(3,1);
rec.cmlAccepted=false;
rec.cmlScoreGain=NaN;
rec.crbMse=NaN;
rec.crbRmse=NaN;

% ---- ST-BA ----
try
    t=tic;
    [th,~]=estimate_stba_k2_runtime(X,cfg);
    rec.runtimeSec(1)=toc(t);
    th=sort(double(th(:).'));
    [rr,mx]=angle_error_k(th,truth.thetaDeg);
    rec.thetaHat(1,:)=th;
    rec.doaMse(1)=rr^2;
    rec.maxErr(1)=mx;
    rec.status(1)="OK";
catch ME
    rec.failed(1)=true;
    rec.status(1)="ERROR:"+string(ME.identifier);
end

% ---- continuous CML ----
try
    t=tic;
    [th,ci]=cml_equalamp_k2_continuous_polish(X,cfg.p);
    rec.runtimeSec(2)=toc(t);
    th=sort(double(th(:).'));
    [rr,mx]=angle_error_k(th,truth.thetaDeg);
    rec.thetaHat(2,:)=th;
    rec.doaMse(2)=rr^2;
    rec.maxErr(2)=mx;
    rec.cmlAccepted=logical(ci.accepted);
    rec.cmlScoreGain=double(ci.scoreGain);
    if ci.accepted
        rec.status(2)="POLISH_ACCEPTED";
    else
        rec.status(2)="GRID_RETAINED";
    end
catch ME
    rec.failed(2)=true;
    rec.status(2)="ERROR:"+string(ME.identifier);
end

% ---- TMR-VAI independent equation-based specialization ----
try
    t=tic;
    [th,~]=tmr_vai_fully_augmented_estimate( ...
        X,cfg.p,2,sourceGrid,musicGrid,'mean');
    rec.runtimeSec(3)=toc(t);
    th=sort(double(th(:).'));
    [rr,mx]=angle_error_k(th,truth.thetaDeg);
    rec.thetaHat(3,:)=th;
    rec.doaMse(3)=rr^2;
    rec.maxErr(3)=mx;
    rec.status(3)="OK";
catch ME
    rec.failed(3)=true;
    rec.status(3)="ERROR:"+string(ME.identifier);
end

% ---- paired exact-model conditional CRB ----
if isfield(truth,'sigma2') && isfinite(truth.sigma2) && truth.sigma2>0
    try
        c=coherent_equalamp_constrained_crb_k2( ...
            cfg,truth.thetaDeg,truth.beta,truth.s,truth.sigma2);
        rec.crbMse=c.rmseDeg^2;
        rec.crbRmse=c.rmseDeg;
    catch
        rec.crbMse=NaN;
        rec.crbRmse=NaN;
    end
end
end
