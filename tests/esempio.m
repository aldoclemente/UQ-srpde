

clc
clear 

addpath(genpath(pwd))

% probabilmente possiamo utilizzare direttamente la funzione:
% lege_eval_multidim(X,k,a,b) che valuta i polinomi di Legendre di k (multiindice (?)) in X punti (restituisce un vettore) su [a,b]^N
%																						Immagino che N dipenda dalle righe di X.

N = 2; rule = @(ii) sum(ii); w = 3; base = 0; % da capire 
Lambda = multiidx_gen(N,rule,w,base);

eval_points = rand(N, 10); % N == 2

% ciclo su tutti i nodi ? 
Psi = zeros(rows(Lambda), size(eval_points)(2));
for i = 1:rows(Lambda)
	Psi(i,:) = lege_eval_multidim(eval_points, Lambda(i,:), 0, 1);
end

%% Problema Giocattolo 1
clc
clear 

addpath(genpath(pwd))
source("../load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 2 * x.^0; u0=0; u1=0;
xx = linspace(0,1,N+1);

A = build_stiffness_1D(xx,mu);
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

%% Problema Giocattolo 2 (ok)
clc
clear all
source("../load.m")
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

%% Problema agli autovalori -u_k'' = k u_k

a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 0 * x.^0; u0=0; u1=0;
A = build_stiffness_1D(a,b,N,mu);
A_ = set_dirichlet_bc_1D(A);

[V, lambda] = eig(A_);

lambda_ex = @(k) -4/h^2*(sin((pi*k)./(2*(N+1)))).^2

[sort(diag(lambda)), h * lambda_ex(1:(N+1))'] %pare ok ...

%% OK adesso cerchiamo di creare la Psi (Psi_{ij} = psi_j ( y_i ) M_h \times q stochastic )
%%									Phi (Phi_{ij} = phi_j ( p_i ) N_h \times n FE )

clc
clear all

q = 2;
n = 25; % giusto per ...

eval_points = rand(1,n);
N = 100;
xx = linspace(0,1,N+1);

Phi = zeros(n, N+1);

for i=1:(N+1)
 Phi(:,i) = lagr_eval(xx(i), xx(setxor(1:numel(xx),[i])), eval_points);
end

size(Phi)


%% check (OK!)
eval_points = xx;
Phi = zeros(N+1, N+1);
for j=1:(N+1)
 Phi(:,j) = lagr_eval(xx(j), xx(setxor(1:numel(xx),[j])), eval_points);
end

%% SR-PDE (locs at nodes)
clc
clear all
source("../load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 0 * x.^0; u0=0; u1=0;
xx = linspace(0,1,N+1);
A = build_stiffness_1D(xx,mu);
M = build_mass_1D(xx);

n = length(xx);
locs =  xx; %rand(n,1);
f_ex = @(x) sin(pi*x);
y = f_ex(locs)' + 0.05^2*rand(n,1);

Phi = zeros(n, N+1);

for j=1:(N+1)
 Phi(:,j) = lagr_eval(xx(j), xx(setxor(1:numel(xx),[j])), locs);
end
lambda = 10.^(-5:0);

for k = 1:length(lambda)
SystemMatrix = [[Phi'*Phi, lambda(k) * A'],
				[lambda(k)*A,		-lambda(k)*M]];

rhs = [Phi'*y; zeros(N+1,1)];

solutions = SystemMatrix\rhs;
f_ = solutions(1:(N+1));

figure() 
plot(xx, f_ex(xx), "color", "red", "linewidth", 2);
grid on
hold on 
plot(xx, f_, "color", "black", "linewidth", 2);

end


%% SR-PDE (random locs)
clc
clear all
source("../load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 0 * x.^0; u0=0; u1=0;
xx = linspace(0,1,N+1);
A = build_stiffness_1D(a,b,N,mu);
M = build_mass_1D(a,b,N);

n = 25;
locs =  rand(n,1);

f_ex = @(x) sin(pi*x);
y = f_ex(locs) + 0.05^2*rand(n,1);

Phi = zeros(n, N+1);

% valutazione è sbagliata per locazioni differenti dai nodi della mesh 
for j=1:(N+1)
 Phi(:,j) = lagr_eval(xx(j), xx(setxor(1:numel(xx),[j])), locs);
end
lambda = 10.^(-5:0);

for k = 1:length(lambda)
SystemMatrix = [[Phi'*Phi, lambda(k) * A'],
				[lambda(k)*A,		-lambda(k)*M]];

rhs = [Phi'*y; zeros(N+1,1)];

solutions = SystemMatrix\rhs;
f_ = solutions(1:(N+1));

figure() 
plot(xx, f_ex(xx), "color", "red", "linewidth", 2);
grid on
hold on 
plot(xx, f_, "color", "black", "linewidth", 2);

end

%% Beh... preoccupiamoci per ora di locazioni spaziali corrispondenti ai nodi della mesh.



 




