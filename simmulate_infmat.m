clear all

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 2; % number of the resonators inside the unit cell
inft = 15; % number of reoccuring unit cells to simmulate an infinite material
N_tot = N*inft; % total number of resonators
spacing = 2; pre_lij = ones(1,N-1).*spacing; %pre_lij(1:2:end) = 1; % spacing between the resonators
lij = repmat([pre_lij,pre_lij(end)],1,inft-1); lij = [lij,pre_lij];
len = 1; pre_li = ones(1,N).*len; % length of the resonator
li = repmat(pre_li,1,inft);


L = sum(pre_li)+sum(pre_lij)+pre_lij(end); % length of the unit cell
xm = [0]; % left boundary points of the resonators
for i = 1:N_tot-1
    xm = [xm,xm(end)+li(i)+lij(i)];
end
xp = xm + li; % right boundary points of the resonators

delta = 0.0001; % small contrast parameter
t = 0; % time
vr = 1; % wave speed inside the resonators
eta_i = 1.5; % perturbance of the wave speed in the i-th resonator
v0 = 1; % wave speed outside the resonators

% implement perturbation
i_pertr = floor(N_tot/2); % indicates which resonator is perturbed
vr = ones(1,N).*vr; 
vr = repmat(vr,1,inft); vr(i_pertr) = vr(i_pertr)+eta_i;


% Settings for modulation
Omega = 0.034; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0.4; % modulation amplitude of kappa
epsilon_rho = 0.6; % modulation amplitude of rho
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end
ks = repmat(ks,inft,1);
rs = repmat(rs,inft,1);
phase_kappa = repmat(phase_kappa,1,inft);
phase_rho = repmat(phase_rho,1,inft);

figure()
hold on
plot(0,0,'b*')
for i = 1:N_tot
    plot([xm(i),xp(i)],zeros(1,2),'r-')
end
for i = 1:inft
    plot(i*L,0,'b*')
end


%% Compute the subwavelength resonant frequencies and eigenmodes for the static case

if N_tot > 1
    C = make_capacitance_finite(N_tot,lij); % capacitance matrix
    [w_cap, v_cap] = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr); % subwavelength resonant frequencies
else
    w_cap = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
end
% w_res = w_cap(real(w_cap)>=0); % positive subwavelength resonant frequencies

fig = figure()
hold on
for i = 1:length(v_cap)
    plot(1:2*N_tot,v_cap(:,i),'-','LineWidth',1, 'color', [.5 .5 .5])
end
legend 
xlabel('$i$','fontsize',18,'Interpreter','latex')
ylabel('$v_i$','fontsize',18,'Interpreter','latex')

% generalised capacitance matrix
C = make_capacitance_finite(N_tot,lij); % capacitance matrix
GCM = delta*diag(vr)^2*diag(1./li)*C;
[V,w_res] = eig(GCM);

% plot eigenvectors
fig = figure()
hold on
for i = 1:length(V)
    plot(1:length(V),V(:,i),'-','DisplayName',['$\omega=$ ', num2str(w_res(i, i))],'LineWidth',1, 'color', [.5 .5 .5])
end
legend 
xlabel('$i$','fontsize',18,'Interpreter','latex')
ylabel('$v_i$','fontsize',18,'Interpreter','latex')
close(fig)


%% Compute the subwavelength resonant frequencies and eigenmodes for the time-modulated case

rhot = @(t) 1./(1+epsilon_rho*cos(Omega*t+phase_rho));
sqrtkappat = @(t) sqrt(1./(1+epsilon_kappa*cos(Omega*t+phase_kappa)));
w3 = @(t) Omega.^2.*epsilon_kappa.*(epsilon_kappa.*sin(Omega*t+phase_kappa).^2+2.*epsilon_kappa.*cos(Omega*t+phase_kappa).^2+2.*cos(Omega*t+phase_kappa))/(2*(epsilon_kappa.*cos(Omega*t+phase_kappa)+1).^(3/2));

C = make_capacitance_finite(N_tot,lij); % capacitance matrix
M = @(t) makeM(t, delta, vr,li, C, rhot, sqrtkappat, w3);
[TOUT, X_bloch,V_mode,V,D] = Hill2BlochModes(M,T,100);

% plot eigenmodes
figure()
hold on
m = 1; % fix time
for n = 1:2*N_tot
    plot(1:2*N_tot,(reshape(X_bloch(m,n,:),[1,size(X_bloch(m,n,:),3)])),'-','LineWidth',1, 'color', [.5 .5 .5])
end







