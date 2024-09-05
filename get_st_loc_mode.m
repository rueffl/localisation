function [W_ts,ts] = get_st_loc_mode(Omega,epsilon_kappa,epsilon_rho,phase_kappa,phase_rho,cs,vr,delta,li,C,t0,k_tr,N,v_init)
%GET_ST_LOC_MODE Computes and plots the space-time localised mode
%   Omega:          frequency of kappa_i and rho_i
%   epsilon_kappa:  modulation amplitude of kappa
%   epsilon_rho:    modulation amplitude of rho
%   phase_kappa:    modulation phase shift of kappa
%   phase_rho:      modulation phase shift of rho
%   cs:             defect coefficients of kappa
%   vr:             wave speed inside the resonators
%   delta:          contrast parameter
%   li:             length of resonators
%   C:              capacitance matrix for fixed alpha
%   t0:             time defect scaling
%   k_tr:           truncation parameter
%   N:              total number of resonators  
%   v_init:         initial guess to solve the ODE

    GCM = delta.*diag(vr)^2*diag(1./li)*C; % generalised capacitance matrix
    invkappat = @(t) diag(ones(1,N)+epsilon_kappa.*cos(Omega.*t+phase_kappa));
    rhot = @(t) diag((ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)).^(-1)); 
    invrhot = @(t) diag(ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)); 
    ct = @(t) diag(cs*exp(-(t/t0).^2));

    dt_invkappat = @(t) diag(-Omega*epsilon_kappa.*sin(Omega.*t+phase_kappa));
    dt_rhot = @(t) diag((Omega*epsilon_rho.*sin(Omega.*t+phase_rho)./((ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)).^2))); 
    dt_ct = @(t) diag(-2.*t*cs*exp(-t.^2)./(t0^2));

    dtdt_rhot = @(t) diag(epsilon_rho.*Omega^2*(2*epsilon_rho.*sin(Omega.*t+phase_rho).^2+epsilon_rho.*cos(Omega.*t+phase_rho).^2 ...
        +cos(Omega.*t+phase_rho))./((1+epsilon_rho.*cos(Omega.*t+phase_rho)).^3));

    A = @(t) invrhot(t)*(dt_invkappat(t)+dt_ct(t))*dt_rhot(t)+invrhot(t)*(invkappat(t)+ct(t))*dtdt_rhot(t)-GCM;
    B = @(t) invrhot(t)*(dt_invkappat(t)+dt_ct(t))*rhot(t)+invrhot(t)*(invkappat(t)+ct(t))*2.*dt_rhot(t);
    D = @(t) invrhot(t)*(invkappat(t)+ct(t))*rhot(t);

    Z = zeros(N,N);
    I = eye(N,N);

    LM = @(t) [A(t), B(t); Z, I]; RM = @(t) [Z, D(t); I, Z];
    M = @(t) inv(RM(t))*LM(t);
    MM = @(t,y) M(t)*y;

    [ts, W_ts] = ode45(MM,[0,k_tr*t0],v_init); 
    

end