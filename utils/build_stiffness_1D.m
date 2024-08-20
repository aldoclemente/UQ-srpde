function [A] = build_stiffness_1D(x,mu) %x nodes
	N = length(x) - 1;
	
	%x = linspace(a,b,N+1);
    h = (x(end)-x(1))/N;

    x_mid = x(1:end-1)+h/2;
    %if(isnumeric(mu))
    %	if(size(mu) == 1 && size(mu) == 1)
    %		a = mu * ones(N+1,1);
    %	else 
    %		a = mu;
    %	endif   
    %else 
    %	a = mu(x_mid);
    %endif
    
    A = 2*diag(ones(1,N+1)) - diag(ones(1,N),1) - diag(ones(1,N),-1);
    A = mu/h * A;
    
end
