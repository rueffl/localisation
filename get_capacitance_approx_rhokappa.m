function [w_cap] = get_capacitance_approx_rhokappa(Omega,epsilon_kappa,epsilon_rho,phase_kappa,phase_rho,vr,delta,li,k_tr,C)
%GET_CAPACITANCE_APPROX_RHOKAPPA Computes the capacitance matrix approximation
%   Omega:          frequency of kappa_i and rho_i
%   epsilon_kappa:  modulation amplitude of kappa
%   epsilon_rho:    modulation amplitude of rho
%   phase_kappa:    modulation phase shift of kappa
%   phase_rho:      modulation phase shift of rho
%   vr:             wave speed inside the resonators
%   delta:          contrast parameter
%   li:             length of resonators
%   k_tr:           truncation parameter
%   C:              capacitance matrix for fixed alpha

%     rhot = @(t) [1./(1+(cos(Omega*t+phir).*epsr))];
%     sqrtkappat = @(t) [1./sqrt(1+cos(Omega*t+phik)*epsk)];
%     w3 = @(t) [Omega^2*(4*cos(Omega*t+phik)*epsk+(3+cos(2*(Omega*t+phik)))*epsk)./(8*(1+cos(Omega*t+phik)*epsk).^2)]; % diagonal entries of matrix W_3
%     T = 2*pi/Omega;
%     N = length(li);
% 
%     matrixM = @(t) Mfunc(t,delta,vr,li,C,rhot,sqrtkappat,w3);
%     [eigvals, V] = hill_exp(T,matrixM,N);
%     w_cap = sort(mink(eigvals,2*N));

    GCM = delta*diag(1./li)*C;

    M = 1; % Number of Fourier coefficients of 1/\kappa
    N = size(GCM,1);
    
    R_mod = zeros(2*M+1,N);
    K_mod = zeros(2*M+1,N);
    for i = 1:N        
        R_mod(:,i) = [epsilon_rho/2*exp(-1i*phase_rho(i)); 1; epsilon_rho/2*exp(1i*phase_rho(i))]; % Fourier coefficients of 1/\rho
        K_mod(:,i) = [epsilon_kappa/2*exp(-1i*phase_kappa(i)); 1; epsilon_kappa/2*exp(1i*phase_kappa(i))]; % Fourier coefficients of 1/\kappa
    end

    ns = -k_tr:k_tr;

    NN = 2*k_tr+1;
    O = diag(ns.'*Omega);
    e = ones(NN,1);
    INN = eye(NN);
    IN = eye(N);
    iK = zeros(NN*N);
    iR = zeros(NN*N);
    R = zeros(NN*N);
    for i = 1:N
        Ki = zeros(NN,NN);
        Ri = zeros(NN,NN);
        for m = -M:M
            Ki = Ki+diag(e(1:NN-abs(m))*K_mod(m+M+1,i),m);
            Ri = Ri+diag(e(1:NN-abs(m))*R_mod(m+M+1,i),m);
        end
        Ii = (i-1)*NN+1:i*NN;
        iK(Ii,Ii) = inv(Ki); %% Fourier coefficients of \kappa
        R(Ii,Ii) = Ri; %% Fourier coefficients of 1/\rho
        iR(Ii,Ii) = inv(Ri); %% Fourier coefficients of \rho        
    end
    iRcR = iR*kron(GCM,INN)*R; 
    
    Z = zeros(NN*N);
    mat = -[kron(IN,O), Z; Z, kron(IN,O)]  - 1i*[Z, iK; -iRcR, Z]; % Kroenecker product to get the RHS matrix

    w_cap = sort(eigs(mat,2*N,'smallestabs'),'ComparisonMethod','real'); % The eigenvalues of "mat" are approximately \omega + n\Omega for |n| < N_fouier. Taking the smallest eigenvalues corresponds to n = 0.

end