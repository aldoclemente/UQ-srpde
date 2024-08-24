% NO forcing, NO BC

clc
clear all

cd "sparse-grids-matlab-kit"
addpath(genpath(pwd))
source("../utils/load.m")

% spatial domain ----------------------------------- 
ab = [0,1];
N_el = 50; 
h = (ab(2)-ab(1))/N_el;
xx = linspace(ab(1),ab(2),N_el+1); 
% --------------------------------------------------

% true random field --------------------------------
randn("seed", 0);
rand("seed",1);

N_true = 5; %2;
mu = 1;
%sigma = [0.25, 0.25]; 
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

% 
mass = build_mass_1D(xx);
I_true = sum(mass * f_true); 

observations = f_true + 0.05*abs(range(f_true))*randn(length(xx),1);

plot(xx, f_true, "color", "black", "linewidth", 4);
hold on
scatter(xx, observations, 40, "r", "filled");
% --------------------------------------------------

% for the definition of the random field -----------
N = 2; 
%N = length(sigma); 
sigma = [0.5, 0.45]; 
% sigma = [0.2, 0.2]; 
% --------------------------------------------------

% observations, locations, mesh_nodes, lambda, mu, sigma, y, PDE_rhs = @(x)0 * x.^0, u0=[], u1=[] 
% observations, mesh_nodes, mesh_nodes, lambda, mu, sigma, y, PDE_rhs, u0, u1
%lambda = 1;
lambda = 10.^(-4:0);

[I, ~] = @(y) QoI(observations, xx, xx, lambda, mu, sigma, y, PDE_rhs, u0, u1); 

% for the sparse grid ------------------------------

knots = @(n) knots_CC(n,-sqrt(3),sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

%% loop over w 

w_vec = 2:6; 
exp_I = zeros(size(w_vec)); 
S_old = [];
Sr_old = []; 
evals_old = []; 
counter = 0; 
start_ = time();
for w = [w_vec] %,10]
    counter = counter+1; 
    S = create_sparse_grid(N,w,knots,lev2knots,rule,S_old);
    Sr = reduce_sparse_grid(S);
    [exp_I(counter),evals] = quadrature_on_sparse_grid(I,S,Sr,evals_old,S_old,Sr_old);
    S_old = S; Sr_old = Sr;  evals_old = evals; 
end
end_ = time();
disp(["computation ended: ", mat2str(round( 100*(end_ - start_))/(100*60)), " mins"]);
err_exp_I = abs(exp_I(1:end-1)-exp_I(end)); 

figure
semilogy(w_vec,err_exp_I,'-*','LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('w','FontSize',18,'Interpreter','latex')
ylabel('Error','FontSize',18,'Interpreter','latex')

err_exp_I = abs(exp_I(1:end)-I_true); 
figure
semilogy(w_vec,err_exp_I,'-*','LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('w','FontSize',18,'Interpreter','latex')
ylabel('Error','FontSize',18,'Interpreter','latex')

%% fix w = 4 

w = 4; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

% computing the mean --------------------------------
%[I, ~] = @(y) QoI(observations, xx, xx, 1, mu, sigma, y, PDE_rhs, u0, u1); 

[exp_I, I_on_Sr] = quadrature_on_sparse_grid(I,Sr); % takes in input the function



% computing the variance ----------------------------
int_I2 = zeros(1, length(lambda));
var_I = int_I2;

for k=1:length(lambda)
	int_I2(k) = quadrature_on_sparse_grid(I_on_Sr(k,:).^2,Sr); % int_{\Gamma} I^2 \rho
	var_I(k) = int_I2(k)-exp_I(k)^2;
end


%% no forcing, no BC -------------------------------------------------------------------------------
clc
clear all

source("../utils/load.m")

% spatial domain ----------------------------------- 
ab = [0,1];
N_el = 50; 
h = (ab(2)-ab(1))/N_el;
xx = linspace(ab(1),ab(2),N_el+1); 
% --------------------------------------------------

% true random field --------------------------------
randn("seed", 0);
rand("seed",1);

N_true = 5; %2;
mu = 1;
%sigma = [0.25, 0.25]; 
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

% 
mass = build_mass_1D(xx);
I_true = sum(mass * f_true); 

observations = f_true + 0.05*abs(range(f_true))*randn(length(xx),1);

plot(xx, f_true, "color", "black", "linewidth", 4);
hold on
scatter(xx, observations, 40, "r", "filled");
% --------------------------------------------------

% for the definition of the random field -----------
N = 2; 
%N = length(sigma); 
sigma = [0.5, 0.45]; 
% sigma = [0.2, 0.2]; 
% --------------------------------------------------

% observations, locations, mesh_nodes, lambda, mu, sigma, y, PDE_rhs = @(x)0 * x.^0, u0=[], u1=[] 
% observations, mesh_nodes, mesh_nodes, lambda, mu, sigma, y, PDE_rhs, u0, u1
%lambda = 1;
lambda = 10.^(-4:0);

% [I, ~] = @(y) QoI(observations, xx, xx, lambda, mu, sigma, y) %, PDE_rhs, u0, u1); 
[I, ~] = @(y) QoI(observations, xx, xx, lambda, mu, sigma, y, PDE_rhs, u0, u1); 

% for the sparse grid ------------------------------

knots = @(n) knots_CC(n,-sqrt(3),sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

%% fix w = 4 
w = 4; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

% computing the mean --------------------------------
%[I, ~] = @(y) QoI(observations, xx, xx, 1, mu, sigma, y, PDE_rhs, u0, u1); 

[exp_I, I_on_Sr] = quadrature_on_sparse_grid(I,Sr); % takes in input the function



% computing the variance ----------------------------
int_I2 = zeros(1, length(lambda));
var_I = int_I2;

for k=1:length(lambda)
	int_I2(k) = quadrature_on_sparse_grid(I_on_Sr(k,:).^2,Sr); % int_{\Gamma} I^2 \rho
	var_I(k) = int_I2(k)-exp_I(k)^2;
end

int_I2
var_I


