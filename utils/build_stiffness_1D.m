function [A] = build_stiffness_1D(x,mu) %x nodes
	N = length(x) - 1;
	
	%x = linspace(a,b,N+1);
    h = (x(end)-x(1))/N;

    x_mid = x(1:end-1)+h/2; 
     
    A = 2*diag(ones(1,N+1)) - diag(ones(1,N),1) - diag(ones(1,N),-1);
    A = mu/h * A;
    %u = [0;A\rhs;0]; % adding boundary terms
end
