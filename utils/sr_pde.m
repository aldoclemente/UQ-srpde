function [f] = sr_pde(observations, locations, mesh_nodes, lambda, 
					  mu = 1, PDE_rhs = @(x)0 * x.^0, u0=[], u1=[])
	
	N = length(mesh_nodes) - 1; 
	A = build_stiffness_1D(mesh_nodes,mu);
	M = build_mass_1D(mesh_nodes);
	u = build_rhs_1D(mesh_nodes, PDE_rhs);
	
	n = length(locations);
	
	Phi = build_Phi(mesh_nodes, locations);
	
	f = zeros(length(mesh_nodes), length(lambda));
	
	for k = 1:length(lambda)
		SystemMatrix = [[Phi'*Phi, lambda(k) * A'],
					   [lambda(k)*A,		-lambda(k)*M]];

		rhs = [Phi'*observations; u];
		if( size(u0)(1) != 0 && size(u1)(1) != 0)
			
			SystemMatrix(1, :) = zeros(1, 2*(N+1));
			SystemMatrix(N+1, :) = zeros(1, 2*(N+1));
			
			SystemMatrix(N+2, :) = zeros(1, 2*(N+1));
			SystemMatrix(2*(N+1), :) = zeros(1, 2*(N+1));
			
			
			SystemMatrix(1,1) = 1; rhs(1) = u0;
			SystemMatrix(N+1,N+1) = 1; rhs(N+1) = u1;
			
			SystemMatrix(N+2,N+2) = 1; rhs(N+2) = u0;
			SystemMatrix(2*(N+1),2*(N+1)) = 1; rhs(2*(N+1)) = u1;
		endif
		
		solutions = SystemMatrix\rhs;
		f(:,k) = solutions(1:(N+1));
	end 
end
