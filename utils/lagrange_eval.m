function [evals] = lagrange_eval(x, i, points) % evals the i-th lagrange basis function, centered at x(i), at points 
	N = length(x) - 1;
	h = (x(end) - x(1))/N;
	evals = zeros(1, length(points));
	
	for j=1:length(points)
		if( abs(points(j)-x(i)) < h )
				evals(j) = 1 - sign((x(i)-points(j)))*(x(i)-points(j))/h;
		end
	end
end
