% clear all
format long

% Settings for the structure
k_tr = 1; % truncation parameters as in remark 3.3
N = 2; % number of the resonators inside the unit cell
inft = 100; % number of reoccuring unit cells to simmulate an infinite material
N_tot = N*inft; % total number of resonators
spacing = 2; 
if N > 1
    pre_lij = ones(1,N-1).*spacing; % spacing between the resonators
    lij = repmat([pre_lij,pre_lij(end)],1,inft-1); lij = [lij,pre_lij];
else
    pre_lij = spacing;
    lij = spacing.*ones(1,inft);
end
len = 1; pre_li = ones(1,N).*len; % length of the resonator
li = repmat(pre_li,1,inft);

if N > 1
    L = sum(pre_li)+sum(pre_lij)+pre_lij(end); % length of the unit cell
else
    L = pre_lij+len; % length of the unit cell
end
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

% prepare plot
fig = figure();
hold on
cm = parula(5);

% settings for modulation
Omega = 0.034; % modulation frequency
T = 2*pi/Omega;
phase_kappa_orig = zeros(1,N); % modulation phases of kappa
phase_rho_orig = zeros(1,N); % modulation phases of rho
for i = 1:(N-1)
    phase_kappa_orig(i+1) = pi/i;
    phase_rho_orig(i+1) = pi/i;
end
eps_k = [0,0.8];
eps_r = [0,0.8];

cc = 1;
for epsilon_kappa = eps_k
    for epsilon_rho = eps_r

        % settings for modulation
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
        phase_kappa = repmat(phase_kappa_orig,1,inft);
        phase_rho = repmat(phase_rho_orig,1,inft);

        % generalised capacitance matrix
        C = make_capacitance_finite(N_tot,lij); % capacitance matrix
        GCM = delta*diag(vr)^2*diag(1./li)*C;
        
        if N_tot > 1
            C = make_capacitance_finite(N_tot,lij); % capacitance matrix
            [w_cap,v_cap] = get_capacitance_approx_rhokappa(Omega,epsilon_kappa,epsilon_rho,phase_kappa,phase_rho,vr,delta,li,k_tr,C); % subwavelength resonant frequencies
        else
            w_cap = get_capacitance_approx_spec_im_N1_1D(epsilon_kappa,Omega,len,delta,vr,v0); % subwavelength resonant frequencies
        end
        w_cap = diag(w_cap); %w_cap = w_cap(real(w_cap)>=0);
        v_cap_rev = [];
        for j = 1:2*N_tot
            for i = 1:N_tot
                v_cap_rev(i,j) = sum(abs(v_cap((i-1)*(2*k_tr+1)+1:i*(2*k_tr+1),j)).^2);
            end
        end
        
        loc_measure = zeros(1,N_tot);
        for i = 1:N_tot
            loc_measure(i) = norm(v_cap_rev(i,:),Inf)/norm(v_cap_rev(i,:),2);
        end
        loc_idx = [];
        for i = 1:N_tot
            if loc_measure(i) > (mean(loc_measure)+0.1)
                loc_idx = [loc_idx,i];
            end
        end

        % create plot
        plot(1:length(loc_measure),loc_measure,'.','LineWidth',2,'MarkerSize',10,'Color', cm(cc,:), ...
            'DisplayName', ['$\varepsilon_{\kappa}=$ ',num2str(epsilon_kappa),', $\varepsilon_{\rho}=$ ',num2str(epsilon_rho)])
        cc = cc+1;

    end
end

% labels for plot
xlabel('$i$','fontsize',18,'Interpreter','latex')
ylabel('$\frac{||v_i^{\alpha}||_{\infty}}{||v_i^{\alpha}||_2}$','fontsize',18,'Interpreter','latex','rotation',0)
xlim([1,length(v_cap_rev(:,end))])
dim = [0.2 0.671111111111111 0.203633567216567 0.128888888888889];
% str = {['$\varepsilon_{\kappa}=$ ',num2str(epsilon_kappa)],['$\varepsilon_{\rho}=$ ',num2str(epsilon_rho)]};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18);
legend('Interpreter','latex','FontSize',14)
