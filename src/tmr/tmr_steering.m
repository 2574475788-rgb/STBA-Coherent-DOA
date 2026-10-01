function A = tmr_steering(p,thetaDeg)
%TMR_STEERING Steering convention used in the TVT paper:
% a(theta) = exp(-j*pi*p*sin(theta)) for d=lambda/2.

p = p(:);
thetaDeg = thetaDeg(:).';
A = exp(-1j*pi*p*sind(thetaDeg));
end
