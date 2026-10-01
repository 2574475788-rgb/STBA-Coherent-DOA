function [rmseTrial,maxErr,errVec] = angle_error_k(thetaHat,thetaTrue)
%ANGLE_ERROR_K Ordered 1-D DOA matching error for equal-cardinality sets.
%
% In one-dimensional angle space with squared-error loss, sorting both sets
% gives the optimal monotone matching.

thetaHat=sort(thetaHat(:).');
thetaTrue=sort(thetaTrue(:).');

assert(numel(thetaHat)==numel(thetaTrue), ...
    'Estimated and true DOA sets must have equal cardinality.');

errVec=thetaHat-thetaTrue;
rmseTrial=sqrt(mean(errVec.^2));
maxErr=max(abs(errVec));
end
