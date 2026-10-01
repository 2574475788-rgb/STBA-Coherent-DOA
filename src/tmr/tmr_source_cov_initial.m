function info = tmr_source_cov_initial(Rx,p,thetaGrid)
%TMR_SOURCE_COV_INITIAL Eq. (17) of Zhang et al.
%
% A_dag(Theta) = A^H (A A^H)^(-1), since A is wide.
% R_s_hat      = A_dag R_x A_dag^H.

Agrid = tmr_steering(p,thetaGrid);
G = Agrid*Agrid';

% Use a linear solve rather than explicit inv; mathematically identical.
Adag = Agrid' / G;

Rs = Adag*Rx*Adag';
Rs = (Rs+Rs')/2;

info.Agrid = Agrid;
info.Adag = Adag;
info.Rs = Rs;
info.rcondGram = rcond(G);
end
