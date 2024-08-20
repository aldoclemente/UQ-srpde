clc
clear 

cd "../sparse-grids-matlab-kit"
addpath(genpath(pwd))

source("../utils/load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 2 * x.^0; u0=0; u1=0;
sigma = [0,0]; y = [0,0];
xx = linspace(0,1,N+1);

A = build_stiffness_rf_1D(xx,mu, sigma, y);
A_ = set_dirichlet_bc_1D(A);
rhs = build_rhs_1D(xx, PDE_rhs);
rhs = [u0, rhs', u1]';
u_ = A_\rhs;

u_ex = @(x) x.*(1-x);
plot(xx, u_ex(xx), "color", "red", "linewidth", 2);
grid on
hold on 
plot(xx, u_, "color", "black", "linewidth", 2);
max(abs(u_ex(xx)' - u_))

