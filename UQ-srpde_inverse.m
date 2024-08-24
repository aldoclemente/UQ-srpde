clc
clear all

cd "sparse-grids-matlab-kit"
addpath(genpath(pwd))

source("../utils/load.m")
% spatial domain ----------------------------------- 
ab = [0,1];
N_el = 50; 
h = (ab(2)-ab(1))/N_el;
x_h = linspace(ab(1),ab(2),N_el+1); 
% --------------------------------------------------

% true random field --------------------------------
randn("seed", 0);
rand("seed",1);

K = 2; %2;
mu = 1;
sigma = [0.5, 0.5];
y_true = -sqrt(3) + (2*sqrt(3)*rand(K,1));

u0 = 0;
u1 = 0;
%PDE_rhs = @(x) pi^2 * sin(pi*x); # ones(size(x));

PDE_rhs = @(x) ones(size(x));

stiff_true = build_stiffness_rf_1D(x_h, mu, sigma, y_true);
rhs_true = build_rhs_1D(x_h, PDE_rhs);
% set dirichlet BC
stiff_true = set_dirichlet_bc_1D(stiff_true);
rhs_true = [u0, rhs_true(2:end-1)', u1]';
% solve 
f_true = stiff_true\rhs_true;

% generate noisy data ------------------------------
sigma_eps = 0.01;
observations = f_true + sigma_eps*randn(rows(f_true),1);

setenv ("OCTAVE_LATEX_DEBUG_FLAG", "1")

figure
scatter(x_h,observations,'*','LineWidth',2)
hold on 
plot(x_h, f_true,'LineWidth',2)
legend("obs",'$f(x,y^*)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

% construct the surrogate model --------------------
lambda = linspace(1e-4,1,5);

knots = @(n) knots_CC(n,-sqrt(3), sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

N = 2;
w = 5; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

f_star = zeros(N_el + 1, length(lambda));
for k=1:length(lambda)
	f = @(y) rf_sr_pde(observations, x_h, x_h, lambda(k), mu, sigma, y, PDE_rhs, u0, u1);  
	f_on_Sr = evaluate_on_sparse_grid(f,Sr); % u_on_Sr has dimension K x Sr.size:
                                         % the ith row contains u(x_i,Sr.knots)
	f_approx = @(y) interpolate_on_sparse_grid(S,Sr,f_on_Sr,y); 

	LS = @(y)  sum( (observations - f_approx(y)).^2 );
	
	%initial guess
	y_start = [0;0]; 

	% minimization of the NLL to find the MAP
	y_MAP = fminsearch(@(y) LS(y),y_start); 
	
	f_star(:,k) = f_approx( y_MAP );
end

line_color = ["y", "g", "m", "c", "b"];
legStr = ["exact"];

figure()
plot(x_h, f_true, "color", "red", 
	 "linestyle", "--","linewidth", 3);
hold on 
grid on
for k=1:length(lambda) 
	hold on 	
	plot(x_h, f_star(:,k), "color", line_color(k), "linewidth", 2);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );


imgdir = "../imgs/";
if(!exist(imgdir))
	mkdir(imgdir)
end


legStr = ["exact"];
figname = [imgdir "img_1.jpg"];
figure()
plot(x_h, f_true, "color", "red", 
	 "linestyle", "--","linewidth", 3);
hold on 
grid on
for k=3:length(lambda) 
	hold on 	
	plot(x_h, f_star(:,k), "color", line_color(k), "linewidth", 3);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );
print(figname)
