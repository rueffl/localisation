function out = makeM(t, delta, vr,li, C, rhot, sqrtkappat, w3)
%MAKE_M  Computes the matrix M in the capacitance matrix approximation ODE-formula
%   Output: returns function handle: coefficient of degree one equation associated to (11.4.10) in Erik's thesis (p. 268)
%   Input:  t               scalar  time variable
%           delta           scalar  High contrast parameter \delta, delta = rho_b/rho_0
%           kr              vector  wave speed inside the resonators 
%           li              scalar  length of resonators
%           C               matrix  capacitance matrix, UNCERTAIN: renormalized or not? 
%           rhot            function handle  (time-dependent) density inside the resonators, vector-valued (NOT matrix-valued)
%           sqrtkappat      function handle  square root of the (time-dependent) bulk modulus inside the resonators, vector-valued (NOT matrix-valued)
%           w3              funciton handle  (time-dependent) vector-valued function (NOT matrix-valued) given by (w3(t))_i = sqrt(kappa_i(t))/2 * d/dt kappa_i'(t)/kappa_i(t)^(3/2), see page 268 in Eriks thesis
%
%   Output: M(t)    matrix  M(t) = delta*kappar/rhor * W1(t) * C * W2(t) + W3(t)
%

    D = diag(li);
    V = diag(vr.^2);
    R = diag(rhot(t));
    Rinv = diag(1./rhot(t));
    K = diag(sqrtkappat(t));
    W3 = diag(w3(t));

    out = delta*V*inv(D)*K*R*C*K*Rinv + W3;

end