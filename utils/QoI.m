function [integral,f] = QoI(observations, locations, mesh_nodes, lambda, 
					  mu, sigma, y, PDE_rhs = @(x)0 * x.^0, u0=[], u1=[])
	
	f = rf_sr_pde(observations, locations, mesh_nodes, lambda, mu, sigma, y, PDE_rhs, u0, u1);
	
	N = length(mesh_nodes)-1; 
    h = (mesh_nodes(end)-mesh_nodes(1))/N;
    
	% trapezoidal rule
	%integral = [h/2, h*ones(1,N-1), h/2]*f; 
	mass = build_mass_1D(mesh_nodes);
	integral = sum( mass * f);
end
