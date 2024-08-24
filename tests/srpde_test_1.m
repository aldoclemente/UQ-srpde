%% SR-PDE locations at nodes (no forcing, no BC)
clc
clear all

cd "../sparse-grids-matlab-kit"
addpath(genpath(pwd))

source("../utils/load.m")
randn("seed", 0);
rand("seed",1);

a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 0*x.^0; u0=[]; u1=[];
xx = linspace(0,1,N+1);
n = length(xx);
locs =  xx; %rand(n,1);
f_ex = @(x) sin(pi*x);
y = f_ex(locs)' + 0.05^2*rand(n,1);

lambda = 10.^(-4:0);

f_ = sr_pde(y, xx, xx, lambda, mu);

line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_ex(xx), "color", "red", 
	 "linestyle", "--","linewidth", 3);
grid on
for k=1:length(lambda) 
	hold on 	
	plot(xx, f_(:,k), "color", line_color(k), "linewidth", 3);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% SR-PDE locations at nodes (NO forcing, + BC)
clc
clear all
source("../utils/load.m")
randn("seed", 0);
rand("seed",1);

a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) 0*x.^0; u0=0; u1=0;
xx = linspace(0,1,N+1);
n = length(xx);
locs =  xx; %rand(n,1);
f_ex = @(x) sin(pi*x);
y = f_ex(locs)' + 0.05^2*rand(n,1);

lambda = 10.^(-4:0);

f_ = sr_pde(y, xx, xx, lambda, mu, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_ex(xx), "color", "red", 
	 "linestyle", "--","linewidth", 3);
grid on
for k=1:length(lambda) 
	hold on 	
	plot(xx, f_(:,k), "color", line_color(k), "linewidth", 3);
	legStr = [legStr; mat2str(lambda(k))];
end
legend( legStr );
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% SR-PDE locations at nodes (+ forcing, NO BC)
clc
clear all
source("../utils/load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) pi^2*sin(pi*x); u0=[]; u1=[];
xx = linspace(0,1,N+1);
n = length(xx);
locs =  xx; %rand(n,1);
f_ex = @(x) sin(pi*x);
y = f_ex(locs)' + 0.05^2*rand(n,1);

lambda = 10.^(-4:0);

f_ = sr_pde(y, xx, xx, lambda, mu, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_ex(xx), "color", "red", 
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

%% SR-PDE locations at nodes (+ forcing, + BC)
clc
clear all
source("../utils/load.m")
a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) pi^2*sin(pi*x); u0=0; u1=0;
xx = linspace(0,1,N+1);
n = length(xx);
locs =  xx;
f_ex = @(x) sin(pi*x);
y = f_ex(locs)' + 0.05^2*rand(n,1);

lambda = 10.^(-4:0);

f_ = sr_pde(y, xx, xx, lambda, mu, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_ex(xx), "color", "red", 
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

%% SR-PDE random locs (+ forcing, + BC)
clc
clear all
source("../utils/load.m")
randn("seed", 0);
rand("seed",1);

a = 0; b=1; N=100; mu=1; PDE_rhs = @(x) pi^2*sin(pi*x); u0=0; u1=0;
xx = linspace(0,1,N+1);
n = 50;
locs =  rand(n,1); 
f_ex = @(x) sin(pi*x);

n = 50;
locs =  rand(n,1);

y = f_ex(locs')' + 0.05*range(f_ex(xx)) * randn(n,1);

figure()
scatter(locs, y, "color", "b")
hold on
plot(xx, f_ex(xx), "linewidth", 3)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

lambda = 10.^(-4:0);

f_ = sr_pde(y, locs, xx, lambda, mu, PDE_rhs, u0, u1);
line_color = ["b", "g", "y", "c", "m"];

legStr = ["exact"];
figure()
plot(xx, f_ex(xx), "color", "red", 
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

