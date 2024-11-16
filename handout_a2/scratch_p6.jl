using LinearAlgebra


function power_method_symmetric(A, tol)

        n = size(A,1);
        v = ones(n);

        lambda = 1;
        r = A*v - lambda*v;

        #for symmetric matrices, the bound evaluates to ||r||_2 /||x||_2
        bound = norm(r)/(norm(v));

        while bound >= tol

                v_tilde = v/norm(v);
                v = A*v_tilde;

                lambda = v[1]/v_tilde[1];

                r = A*v - lambda*v;

                bound = norm(r)/norm(v);


        end

        return lambda, v/norm(v)
end
"""
Compute the eigenpairs of the k extremal eigenvalues (i.e., the k eigenvalues of real
symmetric matrix A with the greatest absolute value).

Inputs:
    A: nxn real symmetric matrix
    k: number of extremal (i.e., maximum magnitude) eigenvalues to compute
    tol: error tolerance; i.e., each eigenvalue estimate
         must be within tol of a true eigenvalue

Outputs:
    λ: vector of k real elements containing the estimates of the extremal eigenvalues;
       λ[i] contains the ith largest eigenvalue by absolute value
    V: nxk matrix where V[:, i] is the eigenvector for the ith largest eigenvalue
       by absolute value

"""


function extremal_eigenpairs(A,k,tol)
	n = size(A,1);
	lambda_vect = Vector{Float64}(undef,k);
	V = Matrix{Float64}(undef,n,k);

	A_mod = copy(A);

	for i=1:k
		lambda, eigvect = power_method_symmetric(A_mod,tol);
		lambda_vect[i] = lambda; 

		V[:,i] = eigvect;
		


		A_mod -= lambda*eigvect*transpose(eigvect);

	end


	return lambda_vect, V;


end

A = [2.0 -1.0  0.0;
    -1.0  2.0 -1.0;
     0.0 -1.0  2.0];

k = 2
tol = 1e-6

# Call your extremal_eigenpairs function
λ, V = extremal_eigenpairs(A, k, tol)

# Print results
println("Computed eigenvalues: ", λ)
println("Computed eigenvectors:")
for i in 1:k
    println(V[:, i])
end


