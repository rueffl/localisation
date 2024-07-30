function out = Mfunc(t,delta,vr,li,C,rhot,sqrtkappat,w3)
% MFUNC     Constructs the matrix M in the capacitance approximation formula
%   t:          time coordinate
%   delta:      contrast parameter
%   vr:         wave speed inside the resonators
%   li:         length of resonators
%   C:          capacitance matrix
%   rhot:       \rho_i(t), anonymous function
%   sqrtkappat: \sqrt{\kappa_i(t)}, anonymous function
%   w3:         diagonal entries of the matrix W_3

    Rho = diag(rhot(t));
    Rinv = diag(1./rhot(t));
    K = diag(sqrtkappat(t));
    W3 = diag(w3(t));
    out = delta/vr^2*inv(diag(li))*K*Rho*C*K*Rinv + W3;

end