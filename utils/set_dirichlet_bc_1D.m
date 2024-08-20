function [A_bc] = set_dirichlet_bc_1D(A)
	A_bc = A;
	A_bc(1,:) = zeros(1, rows(A));
	A_bc(rows(A), :) = zeros(1, rows(A));
	A_bc(1,1) = 1;
	A_bc(end,end) = 1;
end
