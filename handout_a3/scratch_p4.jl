using LinearAlgebra


function subdivide(a, b, ω, tol)

	n = 2;
	
	error_bound = abs(ω)/(4*(n+1))*(b-a)^(n+1);


	while error_bond >= tol

		n += 1;

		error_bound = abs(ω^(n))/(4*(n+1))*((b-a)/n)^(n+1);


	end


    	return n
end

"""
Computes Chebyshev nodes in the interval [a, b] for the function
cos(ω*x) for a maximum absolute error of tol.

Inputs
    a: lower boundary of the interpolation interval
    b: upper boundary of the interpolation interval
    ω: frequency of cos(ω*x)
    tol: maximum absolute error
Output
    x: distinct Chebyshev nodes	 
"""
function chebyshev_nodes(a, b, ω, tol)



	#finding n 
	

	n = 2;

	bound = (b-a)^(n+1)/(2^(2*n+1)) * abs(ω^(n)) * 1/(factorial(n+1));

	while bound >= tol

		n +=1;

		bound = (b-a)^(n+1)/(2^(2*n+1)) * abs(ω^(n)) * 1/(factorial(n+1));

	end



	 h = (b-a)/n;


	 x = Vector{Float64}(undef, n);


	 for i = 1:n

		 x[i] = (1/2)*(a+b) + (1/2)*(b-a)*cos((2*i+1)/(2*n+2)*pi); 


	 end


	 return x
end

# Define parameters
a = 0
b = 1
ω = 2
tol = 1e-6

# Call the function
result = chebyshev_nodes(a, b, ω, tol)

# Print the output
println("Chebyshev nodes: ", result)

