Omega = 0.034;
epsilon_kappa = 0.4;
n = 10;
lambda = 1;

omega_star_p = @(t) (-sqrt(-1).*Omega.*epsilon_kappa.*sin(Omega.*t)-2.*n.*Omega.*(1+epsilon_kappa.*cos(Omega.*t))...
    +sqrt(-1).*sqrt(Omega.^2.*epsilon_kappa.*sin(Omega.*t).^2+4.*lambda.*(1+epsilon_kappa.*cos(Omega.*t))))./(2.*(1+epsilon_kappa.*cos(Omega.*t)));

omega_star_m = @(t) (-sqrt(-1).*Omega.*epsilon_kappa.*sin(Omega.*t)-2.*n.*Omega.*(1+epsilon_kappa.*cos(Omega.*t))...
    -sqrt(-1).*sqrt(Omega.^2.*epsilon_kappa.*sin(Omega.*t).^2+4.*lambda.*(1+epsilon_kappa.*cos(Omega.*t))))./(2.*(1+epsilon_kappa.*cos(Omega.*t)));

ts = linspace(0,1000,1000);
fig = figure()
hold on
plot(ts,omega_star_p(ts),'g-')
plot(ts,omega_star_m(ts),'r-')