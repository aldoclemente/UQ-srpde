%% SR-PDE locations at nodes (+ forcing, + BC)
clc
clear all

cd "../sparse-grids-matlab-kit"
addpath(genpath(pwd))
source("../utils/load.m")

% build true field ------------------------------------------------------------
randn("seed", 0);
rand("seed",1);

a = 0; b=1; N=100;
xx = linspace(0,1,N+1);

N_true = 5; %2;
mu = 1;
sigma = [0.5, 0.45, 0.3, 0.2, 0.1];
y_true = -sqrt(3) + (2*sqrt(3)*rand(N_true,1));

u0 = 0;
u1 = 0;
%PDE_rhs = @(x) pi^2 * sin(pi*x); # ones(size(x));

PDE_rhs = @(x) ones(size(x));

stiff_true = build_stiffness_rf_1D(xx, mu, sigma, y_true);
rhs_true = build_rhs_1D(xx, PDE_rhs);
% set dirichlet BC
stiff_true = set_dirichlet_bc_1D(stiff_true);
rhs_true = [u0, rhs_true(2:end-1)', u1]';
% solve 
f_true = stiff_true\rhs_true;

% -----------------------------------------------------------------------------
locs =  xx; %rand(n,1);
n = size(locs)(1);
observations = f_true + 0.05*range(f_true)*randn(n,1);

lambda = 10.^(-4:0);

N = 2; 
sigma = [0.2, 0.2]; 
y = -sqrt(3) + (2*sqrt(3)*rand(N,1));

f_ = rf_sr_pde(observations, xx, xx, lambda, mu, y, sigma, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_true, "color", "red", 
	 "linestyle", "--","linewidth", 3);
hold on 
grid on
for k=1:length(lambda) 
	hold on 	
	plot(xx, f_(:,k), "color", line_color(k), "linewidth", 2);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

mass = build_mass_1D(xx);
integral = sum( mass * f_ );

%% random locations (+ forcing, + BC)
clc
clear all
source("../utils/load.m")
% build true field ------------------------------------------------------------
randn("seed", 0);
rand("seed",1);

a = 0; b=1; N=100;
xx = linspace(0,1,N+1);

N_true = 5; %2;
mu = 1;
sigma = [0.5, 0.45, 0.3, 0.2, 0.1];
y_true = -sqrt(3) + (2*sqrt(3)*rand(N_true,1));

u0 = 0;
u1 = 0;
%PDE_rhs = @(x) pi^2 * sin(pi*x); # ones(size(x));

PDE_rhs = @(x) ones(size(x));

stiff_true = build_stiffness_rf_1D(xx, mu, sigma, y_true);
rhs_true = build_rhs_1D(xx, PDE_rhs);
% set dirichlet BC
stiff_true = set_dirichlet_bc_1D(stiff_true);
rhs_true = [u0, rhs_true(2:end-1)', u1]';
% solve 
f_true = stiff_true\rhs_true;

% -----------------------------------------------------------------------------
n = 50;
locs =  rand(n,1);

Phi = build_Phi(xx, locs); % Phi_{ij} = phi_j ( p_i );

observations = zeros(n,1);
% f(x) = sum f_i Phi_i(x)
for i=1:rows(Phi)
	observations(i) =  Phi(i,:) * f_true; 
end

observations = observations + 0.05*range(f_true) * randn(n,1);

figure()
scatter(locs, observations, "color", "b")
hold on
plot(xx, f_true, "linewidth", 3)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

lambda = linspace(1e-3, 1, 5);

Q = 5; 
sigma = [0.5, 0.45, 0.3, 0.2, 0.1]; 
%y = -sqrt(3) + (2*sqrt(3)*rand(Q,1));
y = y_true(1:Q) + 0.05 * randn(Q,1);


f_ = rf_sr_pde(observations, locs, xx, lambda, mu, y, sigma, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_true, "color", "red", 
	 "linestyle", "--","linewidth", 3);
hold on 
grid on
for k=3:length(lambda) 
	hold on 	
	plot(xx, f_(:,k), "color", line_color(k), "linewidth", 2);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

mass = build_mass_1D(xx);
integral = sum( mass * f_ );

