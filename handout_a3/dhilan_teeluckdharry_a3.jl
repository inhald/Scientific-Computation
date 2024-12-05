""" 
Computes the coefficients of Newton's interpolating polynomial. 
    Inputs 
        x: vector with distinct elements x[i] 
        y: vector of the same size as x 
    Output 
        c: vector with the coefficients of the polynomial
"""
function newton_int(x, y) 

	#init
        n = length(x);
        c = copy(y);

	#each coefficient is determined by the current
	#and previous indices of c
        for i = 1:n-1

                for j = n:-1:i+1

                        c[j] = (c[j] - c[j-1])/(x[j]-x[j-i]);

                end


        end

	#return 
        return c


end

"""
Evaluates a polynomial with Newton coefficients c 
defined over nodes x using Horner's rule on the points in X.
Inputs 
    c: vector with n coefficients 
    x: vector of n distinct points used to compute c in newton_int 
    X: vector of m points 
Output 
    p: vector of m points
"""
function horner(c, x, X)
	#init dimensions
	n = length(c);
	m = length(X);


	p = Vector{Float64}(undef,m);

	#for each point in X
	for i = 1:m


		#use horner's rule to evaluate
		#the polynomial in linear time
		p_m = c[n];

		for j = n-1:-1:1

			p_m = c[j] + p_m*(X[i]-x[j]); 


		end

		p[i] = p_m;

	end



	#return 
    	return p
end

"""
Computes the number of equally spaced points to use for 
interpolating cos(ω*x) on interval [a, b] for an absolute
error tolerance of tol.

Inputs
    a: lower boundary of the interpolation interval
    b: upper boundary of the interpolation interval
    ω: frequency of cos(ω*x)
    tol: maximum absolute error 
Output
    n: number of equally spaced points to use 	 
"""
function subdivide(a, b, ω, tol)
	
	#init
	n = 2;
	
	#bound for n=1
	h = (b-a)/n;
	error_bound = abs(ω)^(n)/(4*(n+1))*(h)^(n+1);


	#finding n
	while error_bound >= tol

		n+=1;

		h = (b-a)/n;


		error_bound = abs(ω^(n))/(4*(n+1))*(h)^(n+1);


	end


	#for p_n-1
    	return n-1;
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

	#init for n =1
	n = 1;

        bound = (b-a)^(n)/(2^(2*n-1)) * abs(ω^(n)) * 1/(factorial(n));

	#finding n
        while bound >= tol

		n+=1;
                bound = (b-a)^(n)/(2^(2*n-1)) * abs(ω^(n)) * 1/(factorial(n));
		#bound = abs(ω^n)*((b-a)/2)^n/(factorial(n)*2^(n-1));

        end

	#for p_n-1
	n -=1;

         h = (b-a)/n;


         x = Vector{Float64}(undef, n);

	 #pushing nodes
         for i = 1:n

                 x[i] = (1/2)*(a+b) + (1/2)*(b-a)*cos((2*i-1)/(2*n)*pi);


         end


	 #return
         return x
end

"""
Compute the integral ∫f(x)dx over [a, b] with the composite trapezoidal 
rule using r subintervals.

Inputs:
    
    f: function to integrate
    a: lower bound of the definite integral
    b: upper bound of the definite integral
    r: number of subintervals
"""
function composite_trapezoidal_rule(f, a, b, r)

	#init
	h = (b-a)/r; 
	approximate_integral = h/2*(f(a)+f(b));

	#applying composite trapezoidal rule
	for i = 1:r-1
		approximate_integral = approximate_integral + h*f(a+h*i);
	end


	#return
    	return approximate_integral
end

"""
Compute the integral ∫f(x)dx over [a, b] with the composite midpoint 
rule using r subintervals.

Inputs:
    
    f: function to integrate
    a: lower bound of the definite integral
    b: upper bound of the definite integral
    r: number of subintervals
"""
function composite_midpoint_rule(f, a, b, r)
	#init
	h = (b-a)/r; 

	approximate_integral = 0;

	#applying composite midpoint rule
	for i = 1:r
		approximate_integral = approximate_integral + f(a+h*i-h/2);

	end

	approximate_integral = approximate_integral*h;


	#return
    	return approximate_integral
end

"""
Compute the integral ∫f(x)dx over [a, b] with the composite Simpson's 
rule using r subintervals. Note that r must be even because each 
application of Simpson's rule uses a subinterval of length 2*(b-a)/r.
In other words, the midpoints used by the basic Simpson's rule are 
included in the r+1 points on which we evaluate f(x).

Inputs:
    
    f: function to integrate
    a: lower bound of the definite integral
    b: upper bound of the definite integral
    r: even number of subintervals
"""
function composite_simpsons_rule(f, a, b, r)

	#init
	approximate_integral = (f(a)+f(b));
	h = (b-a)/r;


	#sum over even indices
	for i = 2:2:(r-2)
		xi = a+i*h;
		approximate_integral += 2*f(xi);
	end

	#sum over odd indices
	for i= 1:2:(r-1)
		xi= a + i*h
		approximate_integral +=  4*f(xi);
	end

	approximate_integral = approximate_integral * (h/3);
	
	#return
    	return approximate_integral
end

"""
Compute the integral ∫f(x)dx over [a, b] with the adaptive Simpson's 
rule. Return the approximate integral along with the nodes (points) x 
used to compute it.  

Inputs:
    
    f: function to integrate
    a: lower bound of the definite integral
    b: upper bound of the definite integral
    tol: maximum error we can approximately tolerate (i.e., |I_f - Q| <≈ tol)
    max_depth: maximum number of times this function should be recursively called

Returns:
    approximate_integral: the value of the integral ∫f(x)dx over [a, b]
    x: vector containing the nodes which the algorithm used to compute approximate_integral
"""

function recursive_asr(f, a, b, tol, fa, fm, fb, depth, max_depth, x_values)

	#init
    	m = (a + b)/2
    	h = b - a

    	S = h/6 * (fa + 4*fm + fb)

     	#left and right halves
    	left_m = (a + m)/2
    	right_m = (m + b)/2

    	flm = f(left_m)
    	frm = f(right_m)

    	#new nodes
    	push!(x_values, left_m)
    	push!(x_values, right_m)

    	# simpson's rule over [a, m] and [m, b]
    	S_left = h/12 * (fa + 4*flm + fm)
    	S_right = h/12 * (fm + 4*frm + fb)

    	S2 = S_left + S_right

    	# Error estimate
    	if depth >= max_depth || abs(S2 - S)/15 <= tol
        	return S2
    	else
        	# Recurse on [a, m] and [m, b]
        	left_integral = recursive_asr(f, a, m, tol, fa, flm, fm, depth+1, max_depth, x_values)
        	right_integral = recursive_asr(f, m, b, tol, fm, frm, fb, depth+1, max_depth, x_values)
        	# Combine results
        	integral = left_integral + right_integral
        	return integral
    	end
end



function adaptive_simpsons_rule(f, a, b, tol, max_depth)
	
	#init
	fa = f(a)
    	fb = f(b)
    	m = (a + b)/2
    	fm = f(m)
    	depth = 0

    	x_values = Set{Float64}()
    	push!(x_values, a)
    	push!(x_values, m)
    	push!(x_values, b)

    	# call the recursive function
    	integral = recursive_asr(f, a, b, tol, fa, fm, fb, depth, max_depth, x_values)

	#return
    	x = collect(x_values)
    	sort!(x)
    	return integral, x

end