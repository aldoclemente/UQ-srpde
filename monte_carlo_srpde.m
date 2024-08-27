
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

K = 2;
mu = 1;
sigma = [0.5, 0.5];
%y_true = -sqrt(3) + (2*sqrt(3)*rand(K,1));
y_true = [0.9,-1.1];

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
u_true = stiff_true\rhs_true;

% number of observations
n = 80;
locs =  rand(n,1);
Phi = build_Phi(x_h, locs); % Phi_{ij} = phi_j ( p_i );
z = zeros(n,1);
% u(x) = sum u_i Phi_i(x)
for i=1:rows(Phi)
	z(i) =  Phi(i,:) * u_true; 
end
z = z + 0.05*range(u_true)*randn(n,1);

% select "best smoothing parameter"
lambda = 10.^linspace(-3,0,10);
y = y_true(1:2)' + 0.1 * randn(2,1);
u_ = rf_sr_pde(z, locs, x_h, lambda, mu, y, sigma, @(x)0*x.^0, u0, u1);

legStr = [];
figure()
for k=1:length(lambda) 
 	plot(x_h, u_(:,k), "linewidth", 3);
 	hold on
 	legStr = [legStr; mat2str(lambda(k))];
end
grid on
scatter(x_h, u_true, "r", '*');
legend( legStr ,'interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18);
close()

% monte carlo repetitions
M = [50, 100, 250, 500];

%lambda = linspace(1e-3, 1, 5);
lambda = 0.5e-3;
sigma = [0.5, 0.5];
Q = length(sigma); 
%y = -sqrt(3) + (2*sqrt(3)*rand(Q,1));
y = y_true(1:Q) + 0.05 * randn(Q,1);

exp_u = zeros(N+1,length(M));
var_u = zeros(N+1,length(M));

times_ = zeros(1, length(M));

for m=1:length(M)
	u_ = zeros(N+1, M(m));
	bar_ = waitbar(0/M(m));
	start_ = time();
	for k=1:M(m)
		y = y_true(1:Q)' + 0.1 * randn(Q,1);
	
		%u_(:,k) = rf_sr_pde(z, locs, x_h, lambda, mu, y, sigma, PDE_rhs, u0, u1);
		u_(:,k) = rf_sr_pde(z, locs, x_h, lambda, mu, y, sigma, @(x)0*x.^0, u0, u1);
		msg = [mat2str(round(k/M(m) * 100)), "%"];
		waitbar(k/M(m), bar_, msg)
	end
	end_ = time();
	times_(m) = round(100*(end_ - start_))/(100*60); %mins
	exp_u(:,m) = mean(u_,2);
	var_u(:,m) = 1/(M(m)-1) * ( sum(u_.^2,2) - M(m) * exp_u(:,m).^2);
end
close(bar_);

line_color = ["r", "y", "m", "g","b"];
legStr = ["u"];

imgdir = "../imgs/";
if(!exist(imgdir))
	mkdir(imgdir)
end

figname = [imgdir "monte_carlo_data.jpg"];
figure()
plot(x_h, u_true,"color", line_color(1), 'LineWidth',2)
hold on
scatter(locs,z,"b",'*')
legend('u', "data",'interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
print(figname)
close()

figname = [imgdir "monte_carlo_execution_time.jpg"];
figure()
plot(M, times_,"color", line_color(1), 'LineWidth',2)
xlabel("MC samples")
ylabel("execution time [mins]")
grid on
hold on 
scatter(M, times_, "b", "*")
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
print(figname)
close()

legStr = [];
figname = [imgdir "monte_carlo_mean_field.jpg"];
figure()
hold on
grid on
plot(x_h, exp_u(:,k), "color", line_color(end), "linewidth", 3);
legStr = [legStr; mat2str(M(end))];
plot(x_h, u_true, "color", line_color(1), 'LineWidth',2);
legStr = [legStr; "u"];
legend( legStr ,'interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18);
print(figname)


