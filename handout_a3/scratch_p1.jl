using LinearAlgebra







function newton_int(x, y) 
	n = length(x);
	c = copy(y);

	for i = 1:n-1

		for j = n:-1:i+1
			
			c[j] = (c[j] - c[j-1])/(x[j]-x[j-i]);

		end


	end
	
    	return c
end


x = [1.0, 2.0, 3.0, 4.0]
y = [1.0, 4.0, 9.0, 16.0]

coefficients = newton_int(x, y)
println("Computed coefficients: ", coefficients)


