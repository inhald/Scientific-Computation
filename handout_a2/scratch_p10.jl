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
function newton_optimizer(x0, P, d, tol, max_iters)
	


    	return x_trace
end

