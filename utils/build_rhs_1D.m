function [rhs] = build_rhs_1D(x,PDE_rhs) % a,b,N
	N = length(x) - 1;
	% x = linspace(a,b,N+1);
    h = (x(end)-x(1))/N;

    x_mid = x(1:end-1)+h/2; 
    
    f_mid = PDE_rhs(x_mid);
    f_int = f_mid(2:end-1);  
    rhs = h/2*[f_int(1);f_int(1:end-1)'+f_int(2:end)';f_int(end)];

end
