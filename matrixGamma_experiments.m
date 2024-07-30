clear all
format long

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 1; % number of the resonator
spacing = 2; lij = ones(1,N).*spacing;% lij(1:2:end) = 1; % spacing between the resonators
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
vr = ones(1,N).*vr; 
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
epsilon_kappa = 0.6; % modulation amplitude of kappa
epsilon_rho = 0.4; % modulation amplitude of rho
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

all_ws = linspace(-Omega/2,Omega/2,1000)';
all_alphas = linspace(-pi/L,pi/L,1000)';
all_dets = zeros(length(all_ws),length(all_alphas));

for i = 1:length(all_ws)
    for j = 1:length(all_alphas)
        C = 2/spacing*(1-cos(all_alphas(j)*L));
        Gamma = get_matrixGamma(epsilon_kappa,epsilon_rho,all_ws(i),Omega,C,k_tr);
        all_dets(i,j) = svds(Gamma,1,"smallest");
    end
end

% Create Plot
figure()
hold on
s = surf(all_alphas,all_ws,all_dets)
s.EdgeColor = 'none';
xlabel('$\alpha$',Interpreter='latex')
ylabel('$\omega$',Interpreter='latex')

figure()
hold on
for j = 1:length(all_alphas)
    plot(all_ws,all_dets(:,j),'-','DisplayName',strcat('$\alpha=$  ',num2str(all_alphas(j))))
end
% legend on
