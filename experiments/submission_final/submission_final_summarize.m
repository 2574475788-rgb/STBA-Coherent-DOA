function S = submission_final_summarize(doaMse,maxErr,failed)
%SUBMISSION_FINAL_SUMMARIZE Aggregate per-condition metrics.
%
% doaMse/maxErr/failed: Ntrial x Ncondition x Nmethod.

[nTrial,nCond,nMeth]=size(doaMse);
S=struct();
S.nTrial=nTrial;
S.nCond=nCond;
S.nMeth=nMeth;
S.rmse=nan(nCond,nMeth);
S.p025=nan(nCond,nMeth);
S.p05=nan(nCond,nMeth);
S.medianMax=nan(nCond,nMeth);
S.p90Max=nan(nCond,nMeth);
S.p95Max=nan(nCond,nMeth);
S.p99Max=nan(nCond,nMeth);
S.catGt5=nan(nCond,nMeth);
S.failureRate=nan(nCond,nMeth);

for c=1:nCond
    for m=1:nMeth
        ok=~failed(:,c,m) & isfinite(doaMse(:,c,m)) & isfinite(maxErr(:,c,m));
        S.failureRate(c,m)=1-mean(ok);
        if any(ok)
            S.rmse(c,m)=sqrt(mean(doaMse(ok,c,m)));
            x=maxErr(ok,c,m);
            S.p025(c,m)=mean(x<=0.25);
            S.p05(c,m)=mean(x<=0.5);
            S.medianMax(c,m)=median(x);
            S.p90Max(c,m)=prctile(x,90);
            S.p95Max(c,m)=prctile(x,95);
            S.p99Max(c,m)=prctile(x,99);
            S.catGt5(c,m)=mean(x>5);
        end
    end
end
end
