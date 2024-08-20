function [A] = build_stiffness_rf_1D(x, mu, sigma, y)
	N = length(x) - 1;
	
	h = (x(end)-x(1))/N;
	
	x_mid = x(1:end-1)+h/2; 
    
    N_rv = length(y); 
    what_reg = @(x) max(1,ceil(x*N_rv));
    
    RF = zeros(N_rv,1);
    for i = 1:N_rv
        RF(i) = mu + sigma(i)*y(i);
    end 
    a = RF(what_reg(x_mid));
    
    %A = spdiags([-1/h*[a(2:end-1);0],1/h*(a(1:end-1)+a(2:end)),-1/h*[0;a(2:end-1)]],[-1 0 1],N-1,N-1);
	
	A = diag([a(1); (a(1:end-1)+a(2:end)) ; a(end)]) - diag(a,-1) - diag(a,1);  	
	A = 1/h * A;
end
