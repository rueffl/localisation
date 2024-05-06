function [Gamma] = get_matrixGamma(epsilon_kappa,epsilon_rho,w_res,Omega,C)
%UNTITLED6 Summary of this function goes here
%   epsilon_kappa:  modulation amplitude of kappa
%   epsilon_rho:    modulation amplitude of rho
%   w_res:          resonant frequency
%   Omega:          modulation frequency
%   C:              capacitance matrix

    gm2 = zeros(2*k_tr+1,1);
    gm1 = zeros(2*k_tr+1,1);
    g0 = zeros(2*k_tr+1,1);
    gp1 = zeros(2*k_tr+1,1);
    gp2 = zeros(2*k_tr+1,1);

    for n = -k_tr:k_tr
        gm2(n+k_tr+1) = epsilon_kappa*epsilon_rho*0.25*(w_res+(n-2)*Omega)^2+Omega*epsilon_kappa*epsilon_rho*0.25*(w_res+(n-2)*Omega);
        gm1(n+k_tr+1) = C*epsilon_rho*0.5+Omega*epsilon_kappa*0.5*(w_res+(n-1)*Omega)+(epsilon_kappa+epsilon_rho)*0.5*(w_res+(n-1)*Omega)^2;
        g0(n+k_tr+1) = C+(epsilon_rho*epsilon_kappa*0.5+1)*(w_res+n*Omega)^2;
        gp1(n+k_tr+1) = C*epsilon_rho*0.5-Omega*epsilon_kappa*0.5*(w_res+(n+1)*Omega)+(epsilon_kappa+epsilon_rho)*0.5*(w_res+(n+1)*Omega)^2;
        gp2(n+k_tr+1) = epsilon_kappa*epsilon_rho*0.25*(w_res+(n+2)*Omega)^2-Omega*epsilon_kappa*epsilon_rho*0.25*(w_res+(n+2)*Omega);
    end

    Gamma = diag(gp2,2)+diag(gp1,1)+diag(g0)+diag(gm1,-1)+diag(gm2,-2);

end