using LinearAlgebra



"""
Use Newton's method to solve the nonlinear system of equations described in Problems 7-8.
This should work for Euclidean distance measurements in any dimension n.

Inputs:
    x0: initial guess for the position of the receiver in R^n
    P: nxn matrix with known locations of transmitting beacons as columns
    d: vector in R^n where d[i] contains the distance from beacon P[:, i] to x
    tol: Euclidean error tolerance (stop when norm(F(x)) <= tol)
    max_iters: maximum iterations of Newton's method to try

Returns:
    x_trace: Vector{Vector{Float64}} containing each Newton iterate x_k in R^n.

"""
function compute_jacobian(x,p,n)



	J = zeros(n,n);


	for i = 1:n


		for j = 1:n

			J[i,j] = (x[j] - p[j,i])/norm(x .- p[:,i]);

		end

	end

	return J;


end

function compute_fx(x,p,d,n)


	result = zeros(n);

	for i = 1:n

		result[i] = norm(x .- p[:,i]) - d[i];

	end

	return result;

end



function newton(x0, P, d, tol, max_iters)
	
	n = size(P,1);
	x_trace = Vector{Vector{Float64}}();

	k = 0;

	x = copy(x0);

	fx_k = ones(n);


	#Jacobian

	while k <= max_iters && norm(fx_k) >= tol

		J = compute_jacobian(x,P,n);

		fx_k = compute_fx(x,P,d,n);


		s = J \ -fx_k;

		x += s; 

		push!(x_trace, copy(x));

		k+=1;

	end


    	return x_trace
end
# Test case setup for a 3D system
x0 = [0.5, 0.5, 0.5]  # Initial guess (3D)
P = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 2.0 0.0]  # Beacon positions (3 beacons in 3D)
d = [1.0, sqrt(2), 1.0]  # Distances from beacons
tol = 1e-6  # Convergence tolerance
max_iters = 100  # Maximum iterations

# Run Newton's method
x_trace = newton(x0, P, d, tol, max_iters)

# Print results
println("Final estimated position: ", x_trace[end])
println("Trace of positions: ", x_trace)

