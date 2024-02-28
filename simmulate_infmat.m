clear all

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 2; % number of the resonators inside the unit cell
inft = 5; % number of reoccuring unit cells to simmulate an infinite material
N_tot = N*inft; % total number of resonators
spacing = 2; pre_lij = ones(1,N-1).*spacing; pre_lij(1:2:end) = 1; % spacing between the resonators
lij = repmat([pre_lij(end)/2,pre_lij],1,inft-1); lij = [lij,pre_lij];
len = 1; pre_li = ones(1,N).*len; % length of the resonator
li = repmat(pre_li,1,inft);
L = sum(li)+sum(lij); % length of the unit cell
xm = [lij(end)/2]; % left boundary points of the resonators
for i = 2:N_tot
    xm = [xm,xm(end)+li(i-1)+lij(i-1)];
end
xp = xm + li; % right boundary points of the resonators

delta = 0.0001; % small contrast parameter
t = 0; % time
vr = 1; % wave speed inside the resonators
eta_i = 0.2; % perturbance of the wave speed in the i-th resonator
v0 = 1; % wave speed outside the resonators

% implement perturbation
i_pertr = 2; % indicates which resonator is perturbed
vr = ones(1,N).*vr; vr(i_pertr) = vr(i_pertr)+eta_i;
vr = repmat(vr,1,inft);


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
epsilon_rho = 0.2; % modulation amplitude of rho
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
phase_kappa = repmat(phase_kappa,inft,1);
phase_rho = repmat(phase_rho,inft,1);

% find initial guess of the subwavelength quasifrequencies
if N_tot > 1
    C = make_capacitance_finite(N_tot,lij); % capacitance matrix
    w_cap = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr); % subwavelength resonant frequencies
else
    w_cap = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
end
w_res = w_cap(real(w_cap)>=0); % positive subwavelength resonant frequencies



