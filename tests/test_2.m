clc
clear all

cd "../sparse-grids-matlab-kit"
addpath(genpath(pwd))

source("../utils/load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) (1+pi^2) * sin(pi*x); u0=0; u1=0;

xx = linspace(0,1,N+1);
A = build_stiffness_1D(xx,mu);
M = build_mass_1D(xx,N);
A_ = set_dirichlet_bc_1D(A+M);

rhs = build_rhs_1D(xx, PDE_rhs);
rhs = [u0, rhs', u1]';
u_ = A_\rhs;

u_ex = @(x) sin(pi*x);
figure()
plot(xx, u_ex(xx), "color", "red", "linewidth", 2);
grid on
hold on 
plot(xx, u_, "color", "black", "linewidth", 2);
max(abs(u_ex(xx)' - u_))

