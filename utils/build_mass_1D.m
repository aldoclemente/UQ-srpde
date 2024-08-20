function [M] = build_mass_1D(x) %x nodes
	N = length(x) - 1;
	%x = linspace(a,b,N+1);
    h = (x(end)-x(1))/N;
    
    M = 4*diag(ones(1,N+1)) + diag(ones(1,N),1) + diag(ones(1,N),-1);
    M = h/6 * M;
end
