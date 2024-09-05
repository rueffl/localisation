% clear all
format long

% Settings for the structure
k_tr = 4; % truncation parameters as in remark 3.3
N = 3; % number of the resonators inside the unit cell
inft = 25; % number of reoccuring unit cells to simmulate an infinite material
N_tot = N*inft; % total number of resonators
spacing = 1; 
if N > 1
    pre_lij = ones(1,N-1).*spacing; % spacing between the resonators
    pre_lij(2) = 2;
    lij = repmat([pre_lij,pre_lij(end)],1,inft-1); lij = [lij,pre_lij,pre_lij(end)];
else
    pre_lij = spacing;
    lij = spacing.*ones(1,inft);
end
len = 1; pre_li = ones(1,N).*len; %pre_li(2) = 2;% length of the resonator
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
eta_i = 2; % perturbance of the wave speed in the i-th resonator
v0 = 1; % wave speed outside the resonators

% implement perturbation
i_pertr = floor(N_tot/2); % indicates which resonator is perturbed
vr = ones(1,N).*vr; 
vr = repmat(vr,1,inft); vr(i_pertr) = vr(i_pertr)*eta_i;


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
phase_kappa = repmat(phase_kappa,1,inft);
phase_rho = repmat(phase_rho,1,inft);
cs = 0.*ones(1,N_tot); cs(1,2:N:end) = 20; % cs(1,floor(N_tot/2)) = 200; % temporal defect
t0 = T;

alpha = 0.08;
C = make_capacitance(N_tot,lij,alpha,L);
[w_out,v_pre_out] = get_capacitance_approx(Omega,epsilon_kappa,epsilon_rho,phase_kappa,phase_rho,cs,vr,delta,li,C,t0,k_tr);

v_out = v_pre_out(1:N,:);

figure()
hold on
plot(1:2*N_tot,abs(w_out),'*')

% figure(1)
% hold on
% % subplot(2,2,4)
% plot(real(w_out),imag(w_out),'*')
% xlabel('Re($\omega$)',Interpreter='latex')
% ylabel('Im($\omega$)',Interpreter='latex')

loc_measure = zeros(1,2*N_tot);
for i = 1:2*N_tot
    loc_measure(i) = norm(v_out(:,i),Inf)/norm(v_out(:,i),2);
end
% loc_idx = [];
% for i = 1:N_tot
%     if loc_measure(i) > (mean(loc_measure)+0.1)
%         loc_idx = [loc_idx,i];
%     end
% end

figure()
hold on
plot(1:length(loc_measure),loc_measure,'.','LineWidth',2,'MarkerSize',10)

%% Search for the localised solutions

figure()

for i = 1:2*N_tot

    [W_ts_forward,ts_forward] = get_st_loc_mode(Omega,epsilon_kappa,epsilon_rho,phase_kappa, ...
        phase_rho,cs,vr,delta,li,C,t0,k_tr,N_tot,v_pre_out(:,i));
    
    [W_ts_backward,ts_backward] = get_st_loc_mode(Omega,epsilon_kappa,epsilon_rho,phase_kappa, ...
        phase_rho,cs,vr,delta,li,C,t0,-k_tr,N_tot,v_pre_out(:,i));
    
    ts = [flip(ts_backward);ts_forward];
    num_tsteps = length(ts);
    W_ts = [flip(W_ts_backward);W_ts_forward]; W_ts = W_ts(:,1:N_tot);

    loc_measure_ts = zeros(1,num_tsteps);
    for idx = 1:num_tsteps
        loc_measure_ts(idx) = norm(W_ts(idx,:),Inf)/norm(W_ts(idx,:),2);
    end

    plot(ts,loc_measure_ts,'-')
    ylim([0,1])
    title(strcat('$i=',num2str(i),'$'),Interpreter='latex')
    pause(0.5)

end

%% Compute the localised sol

[W_ts_forward,ts_forward] = get_st_loc_mode(Omega,epsilon_kappa,epsilon_rho,phase_kappa, ...
    phase_rho,cs,vr,delta,li,C,t0,k_tr,N_tot,v_pre_out(:,34));

[W_ts_backward,ts_backward] = get_st_loc_mode(Omega,epsilon_kappa,epsilon_rho,phase_kappa, ...
    phase_rho,cs,vr,delta,li,C,t0,-k_tr,N_tot,v_pre_out(:,34));

ts = [flip(ts_backward);ts_forward];
W_ts = [flip(W_ts_backward);W_ts_forward]; W_ts = W_ts(:,1:N_tot);

% figure()
% ylim([-3.6,4.53])
num_tsteps = length(ts);
x_len = length(W_ts(1,:));
loc_measure_ts = zeros(1,num_tsteps);
for idx = 1:num_tsteps
%     plot(1:x_len,W_ts(idx,:),'-')
% %     ylim([min(min(real(W_ts))),max(max(real(W_ts)))])
%     title(strcat('$t=',num2str(ts(idx)),'$ s'),Interpreter='latex')
%     pause(0.5)
    loc_measure_ts(idx) = norm(W_ts(idx,:),Inf)/norm(W_ts(idx,:),2);
end

figure()
plot(ts,loc_measure_ts,'-')

figure()
c = floor(num_tsteps/15);
for plt_idx = 1:16
    subplot(4,4,plt_idx)
    plot(1:x_len,W_ts((plt_idx-1)*c+1,:),'-')
    ylim([-3,3])
    title(strcat('$t=',num2str(ts((plt_idx-1)*c+1)),'$ s'),Interpreter='latex')
end

