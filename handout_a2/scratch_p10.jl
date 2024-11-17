using LinearAlgebra

"""

Use Newton's method to solve the nonlinear optimization problem described in Problems 9-10.
This should work for Euclidean distance measurements in any dimension n, and any number
    of noisy measurements m.

Inputs:
    x0: initial guess for the position of the receiver in R^n
    P: nxm matrix with known locations of transmitting beacons as columns
    d: vector in R^m where d[i] contains the noisy distance from beacon P[:, i] to x
    tol: gradient Euclidean error tolerance (stop when norm(∇f(x)) <= tol)
    max_iters: maximum iterations of Newton's method to try

Returns:
    x_trace: Vector{Vector{Float64}} containing each Newton iterate x_k in R^n.

"""
function compute_gradf(x, P, d)

	m = size(P)[2];
	n = size(P)[1];

	grad_f = zeros(n);


	for i = 1:m
		diff = x .- P[:,i];
		dist = norm(diff);
		f_i = dist - d[i];

		grad_f .+= 2 * f_i * (diff/dist);


	end


	return grad_f;

end


function compute_hessian(x,P,d)

	m = size(P)[2];
	n = size(P)[1];

	grad_gradf = zeros(n,n);

	

	for i = 1:m

		diff = x .- P[:,i];
		dist = norm(diff);
		f_i = dist - d[i];

		grad_gradf .+=  2*((diff * diff')/dist^2  + f_i * (I/dist  - (diff*diff')/(dist)^3) );

	end


	return grad_gradf; 



end



function newton_optimizer(x0, P, d, tol, max_iters)

	m = size(P)[2];
	n = size(P)[1];

	x_trace = Vector{Vector{Float64}}();
	
	k = 0;


	x = copy(x0);

	grad_f = ones(n);

	while k <= max_iters && norm(grad_f) >= tol
		
		grad_f = compute_gradf(x,P,d);

		grad_gradf = compute_hessian(x,P,d);


		s = grad_gradf \ -grad_f;

		x += s;

		push!(x_trace, copy(x));


		k +=1;

	end
	
    	return x_trace
end

x0 = [0.25, 0.25]                       # Initial guess
P = [0.0 2.0; 0.0 2.0]                # Transmitter positions
d = [1.5, 1.5]                        # Measured distances
tol = 1e-6                            # Tolerance
max_iters = 100                       # Max iterations

x_trace = newton_optimizer(x0, P, d, tol, max_iters)
println(x_trace[end])      
