clear all

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
n = 3; % number of the resonators inside the unit cell
inft = 5; % number of reoccuring unit cells to simmulate an infinite material
N = n*inft; % total number of resonators
spacing = 2; lij = ones(1,N).*spacing; lij(1:2:end) = 1; % spacing between the resonators
perturbation = 0.2; lij(2) = lij(2)-perturbation; lij(3) = lij(3)+perturbation; % add a spatial perturbation
len = 1; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
xm = [lij(end)/2]; % left boundary points of the resonators
for i = 2:N
    xm = [xm,xm(end)+li(i-1)+lij(i-1)];
end
xp = xm + li; % right boundary points of the resonators

delta = 0.0001; % small contrast parameter
t = 0; % time
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

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







