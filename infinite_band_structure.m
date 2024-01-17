% Settings for the structure
k_tr = 2; % truncation parameters as in remark 3.3
N = 2; % number of the resonator
spacing = 1; lij = ones(1,N).*spacing; lij(1:2:end) = 2; % spacing between the resonators
len = 1; li = ones(1,N).*len; % length of the resonator
L = sum(li)+sum(lij); % length of the unit cell
xm = [0,li(1:end-1)+lij(1:end-1)]; % left boundary points of the resonators
xp = xm + li; % right boundary points of the resonators
delta = 0.0001; % small contrast parameter

t = 0; % time
vr = 1; % wave speed inside the resonators
v0 = 1; % wave speed outside the resonators

% Settings for modulation
Omega = 0.03; % modulation frequency
T = 2*pi/Omega;
phase_kappa = zeros(1,N); % modulation phases of kappa
phase_rho = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa(i+1) = pi/i;
    phase_rho(i+1) = pi/i;
end
epsilon_kappa = 0; % modulation amplitude of kappa
epsilon_rho = 0.5; % modulation amplitude of rho
rs = []; % Fourier coefficients of 1/rho
ks = []; % Fourier coefficients of 1/kappa
for j = 1:N
    rs_j = [epsilon_rho*exp(-1i*phase_rho(j))./2,1,epsilon_rho*exp(1i*phase_rho(j))./2];
    ks_j = [epsilon_kappa*exp(-1i*phase_kappa(j))./2,1,epsilon_kappa*exp(1i*phase_kappa(j))./2];
    ks = [ks; ks_j];
    rs = [rs; rs_j];
end

%% Test the function whose roots must be determined

p = @(w,alpha) minev(getMatcalA(N,lij,xm,xp,k_tr,w,Omega,rs,ks,vr,delta,v0,alpha,L));
pts_alpha = 205; pts_w = 205;
alphas = linspace(-pi/L,pi/L,pts_alpha);
ws = linspace(-1*Omega/2,1*Omega/2,pts_w);
ps = zeros(pts_alpha,pts_w);

for ia = 1:pts_alpha
    for iw = 1:pts_w
        ps(iw,ia) = p(ws(iw),alphas(ia));
    end
end

% create plot
figure
s = surf(alphas,ws,real(ps));
s.EdgeColor = 'none';
colorbar

%% Create Band Functions

% Band functions computation
sample_points = 80;
alphas = linspace(-pi/L,pi/L,sample_points);
w_muller = zeros(N,sample_points);


for j = 1:sample_points

    alpha = alphas(j); % quasi periodicity

    % find initial guess of the subwavelength quasifrequencies
    if N > 1
        C = make_capacitance(N,lij,alpha,L); % capacitance matrix
        w_cap = get_capacitance_approx_spec(epsilon_kappa,phase_kappa,Omega,delta,li,v0,vr,C,k_tr); % subwavelength resonant frequencies
    else
        w_cap = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0);
    end
    w_res = w_cap(real(w_cap)>=0); % positive subwavelength resonant frequencies

    % apply Mullers method
    for i = 1:N
        w_muller(i,j) = muller(w_res(i),alpha,N,lij,L,xm,xp,k_tr,Omega,rs,ks,vr,delta,v0);
        while w_muller(i,j) > Omega/2
            w_muller(i,j) = w_muller(i,j)-Omega;
        end
        while w_muller(i,j) < -Omega/2
            w_muller(i,j) = w_muller(i,j)+Omega;
        end
        if real(w_muller(i,j)) < 0
            w_muller(i,j) = -w_muller(i,j);
        end
    end
    [vct,idx] = sort(real(w_muller(:,j)));
    w_muller(:,j) = w_muller(idx,j);

end

% create plot
fig = figure()
hold on
for i = 1:N
    plot(alphas,real(w_muller(i,:)),'-',markersize=8,linewidth=2)
end
xlabel('$\alpha$',fontsize=18,interpreter='latex')
ylabel('$\omega_i$',fontsize=18,interpreter='latex')

