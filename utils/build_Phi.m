function [Phi] = build_Phi(mesh_nodes, locations)
	n = length(locations);
	N = length(mesh_nodes) - 1;
	
	Phi = zeros(n, N+1);
	for j=1:(N+1)
 		%Phi(:,j) = lagr_eval(mesh_nodes(j), mesh_nodes(setxor(1:numel(mesh_nodes),[j])), locations);
 		Phi(:,j) = lagrange_eval(mesh_nodes, j, locations);
	end
end
