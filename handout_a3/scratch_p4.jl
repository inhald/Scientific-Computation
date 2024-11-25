using LinearAlgebra


function subdivide(a, b, ω, tol)

	n = 1;
	
	error_bound = ω/(4*(n+1))*(b-a)^(n+1);


	while error_bond >= tol

		n += 1;

		error_bound = ω^(n-1)/(4*(n+1))*((b-a)/n)^(n+1);


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
	


	

	



    	return x
end


