function [w_out,v_out] = get_capacitance_approx(Omega,epsilon_kappa,epsilon_rho,phase_kappa,phase_rho,cs,vr,delta,li,C)
%GET_CAPACITANCE_APPROX_RHOKAPPA Computes the capacitance matrix approximation
%   Omega:          frequency of kappa_i and rho_i
%   epsilon_kappa:  modulation amplitude of kappa
%   epsilon_rho:    modulation amplitude of rho
%   phase_kappa:    modulation phase shift of kappa
%   phase_rho:      modulation phase shift of rho
%   cs:             defect coefficients of kappa
%   vr:             wave speed inside the resonators
%   delta:          contrast parameter
%   li:             length of resonators
%   k_tr:           truncation parameter
%   C:              capacitance matrix for fixed alpha

    N = length(phase_kappa);
    T = 2*pi/Omega;

    GCM = delta.*diag(vr)^2*diag(1./li)*C; % generalised capacitance matrix

    kappat = @(t) diag((ones(1,N)+epsilon_kappa.*cos(Omega.*t+phase_kappa)).^(-1));
    invkappat = @(t) diag(ones(1,N)+epsilon_kappa.*cos(Omega.*t+phase_kappa));
    rhot = @(t) diag((ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)).^(-1)); 
    invrhot = @(t) diag(ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)); 
    ct = @(t) diag(cs*exp(-t.^2));

    dt_invkappat = @(t) diag(-Omega*epsilon_kappa.*sin(Omega.*t+phase_kappa));
    dt_rhot = @(t) diag((Omega*epsilon_rho.*sin(Omega.*t+phase_rho)./((ones(1,N)+epsilon_rho.*cos(Omega.*t+phase_rho)).^2))); 
    dt_ct = @(t) diag(-2.*t*cs*exp(-t.^2));

    dtdt_rhot = @(t) diag(epsilon_rho.*Omega^2*(2*epsilon_rho.*sin(Omega.*t+phase_rho).^2+epsilon_rho.*cos(Omega.*t+phase_rho).^2 ...
        +cos(Omega.*t+phase_rho))./((1+epsilon_rho.*cos(Omega.*t+phase_rho)).^3));

    A = @(t) invrhot(t)*(dt_invkappat(t)+dt_ct(t))*dt_rhot(t)+invrhot(t)*(invkappat(t)+ct(t))*dtdt_rhot(t)-GCM;
    B = @(t) invrhot(t)*(dt_invkappat(t)+dt_ct(t))*rhot(t)+invrhot(t)*(invkappat(t)+ct(t))*2.*dt_rhot(t);
    D = @(t) invrhot(t)*(invkappat(t)+ct(t))*rhot(t);

    Z = zeros(N,N);
    I = eye(N,N);
    II = eye(2*N,2*N);

    LM = @(t) [A(t), B(t); Z, I]; RM = @(t) [Z, D(t); I, Z];
    M = @(t) inv(RM(t))*LM(t);
    MM = @(t,y) M(t)*y;

    for j = 1:2*N
        [~, wj] = ode45(MM,[0,T],II(:,j)); 
        W(:,j) = wj(end,:)';
%         U = zeros(size(wj));
%         for i = 1:2*N
%             U(:,i) = wj(:,i).*exp()
%         end
    end
    [U, D, V] = eig(W);
    w_out = (log(diag(D))/(1i*T));
    [out_real,ind] = sort(real(w_out),'descend');
    w_out = w_out(ind);
    v_out = U(:,ind);

end


% 
% 
% function [w_out] = get_capacitance_approx(epsilon_kappa,epsilon_rho,li,Omega,phase_rho,phase_kappa,delta,C)
%     T = 2*pi/Omega;
%     N = length(phase_kappa);
%     sqrtkappat = @(t) 1./sqrt(1+epsilon_kappa*cos(Omega*t+phase_kappa));
%     w3 = @(t) Omega^2/4*(1+((epsilon_kappa^2-1)./(1+epsilon_kappa*cos(Omega*t+phase_kappa)).^2));
%     rhot =  @(t) 1./(1 + epsilon_rho*cos(Omega*t+phase_rho));
%     M = @(t) makeM(t,delta,li,C,rhot,sqrtkappat,w3);
%     [w_out, cnd] = hill_exp(T,M,N);
%     [w_out_real,order] = sort(real(w_out),'descend');
%     w_out_imag = imag(w_out(order));
%     w_out = w_out_real + sqrt(-1).*w_out_imag;
% end
% 
% function out = makeM(t,delta,li,C,rhot,sqrtkappat,w3)
%     Rho = diag(rhot(t));
%     Rinv = diag(1./rhot(t));
%     K = diag(sqrtkappat(t));
%     W3 = diag(w3(t));
%     out = delta*diag(1./li)*K*Rho*C*K*Rinv + W3;
% end
% 
% function [w_out, cnd] = hill_exp(T,M,N)
%     W = zeros(2*N,2*N);
%     Z = zeros(N,N);
%     I = eye(N,N);
%     II = eye(2*N,2*N);
%     MM = @(t,y) [Z, I; -M(t), Z]*y;
% %     for j = 1:N
% %         wI0 = zeros(2*N,1); wI0(j) = 1;
% %         wII0 = zeros(2*N,1); wII0(N+j) = 1;
% %         [~, wI] = ode45(MM, [0,T], wI0);
% %         [~, wII] = ode45(MM, [0,T], wII0);
% %         W(:,j) = wI(end,:).';
% %         W(:,N+j) = wII(end,:).';
% %     end
%     for j = 1:2*N
%         [~, wj] = ode45(MM,[0,T],II(:,j)); 
%         W(:,j) = wj(end,:);
%     end
%     [U, D, V] = eig(W);
%     w_out = (log(diag(D))/(1i*T));
%     [out_real,ind] = sort(real(w_out),'descend');
%     w_out = w_out(ind);
%     cnd = cond(U);
% end