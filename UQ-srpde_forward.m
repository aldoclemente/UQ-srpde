
%%% ------------------------------------------------------------------------------------------------
cd "sparse-grids-matlab-kit"
addpath(genpath(pwd))

clc; clear all;
source("../utils/load.m")
% spatial domain ----------------------------------- 
ab = [0,1];
N = 100; 
h = (ab(2)-ab(1))/N;
x_h = linspace(ab(1),ab(2),N+1); 
% --------------------------------------------------

% true random field --------------------------------
randn("seed", 0);
rand("seed",1);

mu = 1;
sigma = [0.5, 0.1];

u0 = 0;
u1 = 0;

PDE_rhs = @(x) ones(size(x));

stiff_true = @(y) build_stiffness_rf_1D(x_h, mu, sigma, y);
rhs_true = build_rhs_1D(x_h, PDE_rhs);

knots = @(n) knots_CC(n,-sqrt(3),sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

% locations
n = 20;
locs = rand(n,1); 
Q = 2;
Phi = build_Phi(x_h, locs); % Phi_{ij} = phi_j ( p_i );

w = 2; 
S = create_sparse_grid(Q,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

obs = zeros(n*Sr.size, 1);
for(i = 1:Sr.size)
	y_i = Sr.knots(:,i) + 0.01*randn(Q,1);
	stiff_i = stiff_true(y_i);
	stiff_i = set_dirichlet_bc_1D(stiff_i);
	u_i = stiff_i\rhs_true;
	obs((1+n*(i-1)):(i*n)) = Phi*u_i+ 0.05*randn(n,1);
end

imgdir = "../imgs/";
if(!exist(imgdir))
	mkdir(imgdir)
end

figname = [imgdir "forward_UQ_data.jpg"];
legStr = [];
figure()
hold on
for(i = 1:Sr.size)
	scatter(locs,obs((1+n*(i-1)):(i*n)),'*','LineWidth',2)
	legStr = [legStr; mat2str(i)];
end
legend(legStr,'interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
print(figname)
close()

u = @(y) rf_sr_pde(obs, repmat(locs,Sr.size,1), x_h, 1, mu, sigma, y, PDE_rhs, u0, u1); 

mass = build_mass_1D(x_h);
QoI = @(y) sum(mass * u(y));
[exp_I, I_on_Sr] = quadrature_on_sparse_grid(QoI, Sr);
int_I2 = quadrature_on_sparse_grid(I_on_Sr.^2, Sr);
var_I = int_I2 - exp_I^2;

% ---
domain = [-sqrt(3), -sqrt(3); sqrt(3), sqrt(3)];
figname = [imgdir, "forward_UQ_sparse_grid_surface.png"];
plot_sparse_grids_interpolant(S, Sr, domain, I_on_Sr);
print(figname)
close()

M = 5000;
y_samples = -sqrt(3) + 2*sqrt(3)*rand(Q,M);

I_vals_rand = interpolate_on_sparse_grid(S, Sr, I_on_Sr, y_samples);
figname = [imgdir, "forward_UQ_hist.png"];
hist(I_vals_rand, nbins=25, 1)
print(figname)
close()
