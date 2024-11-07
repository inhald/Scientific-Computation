using LinearAlgebra

"""
Solve Ax=b iteratively using the Jacobi method with initial iteration x0. 

Inputs:
    A: nxn real matrix 
    b: n-element real vector
    x0: initial iteration
    x_true: ground-truth value used to compute output errs 
    k_max: maximum number of iterations
    res_tol: residual error tolerance criterion; terminate once norm(A*x - b) <= res_tol
Outputs:
    x: n-element real vector containing the solution returned by the Jacobi method
    errs: vector of k+1 Euclidean norm error values, where k is the number of 
          iterations performed. The i-th element errs[i] contains the absolute error 
          in x after i-1 iterations of the Jacobi method have been performed (e.g., 
          errs[1] contains norm(x0 - x_true)).
    bounds: upper bound on errs computed using the convergence theory discussed in lecture. 
"""
function jacobi_method(A, b, x0, k_max)

	n = size(A,1);
	x = Vector{Float64}(undef,n);
	x_prev = copy(x0);
	
	for k = 1:k_max
		for i = 1:n
			sum_aij_xj = 0; 

			for j = 1:n
				if j !== i
					sum_aij_xj += A[i,j]*x_prev[j];
				end
			end

			x[i] = (b[i] - sum_aij_xj)/A[i,i];

		end

		x_prev = copy(x);

	end


    	return x;
end




function ldl_solve(L,d,b)

	#let L^T x = r =>  Dr = s => Ls = b
	#applying forward substitution to find b
	
	n = size(L,1);
	s = Vector{Float64}(undef, size(L,1));
	r = Vector{Float64}(undef, size(L,1));
	x = Vector{Float64}(undef, size(L,1));

	s[1] = b[1];


	for i=2:n
		elem_sum = 0;
		for j = 1:i-1
			elem_sum += L[i,j]*s[j];
		end

		s[i] = b[i] - elem_sum; 
			
	end

	#dividing each element by diagonal element to get r
	
	for i = 1:n
		r[i] = s[i]/d[i];
	end

	#applying back substitution to get x
	
	x[n] = r[n];

	for i = n-1:-1:1

		row_sum = 0; 

		for j = n:-1:i+1
			
			row_sum+= L[j,i]*x[j]
		end

		x[i] = r[i] - row_sum;


	end

	return x; 

end

function ldl_decomposition(A)

	n = size(A,1);
	d = Vector{Float64}(undef,n);
	L = Matrix{Float64}(I,n,n);
	
	d[1] = A[1,1];

	for j = 1:n-1

		for i = j+1:n

			sum_ldl = 0;

			for k = 1:j-1
				sum_ldl += L[i,k]*d[k]*L[j,k];
			end

			L[i,j] = (A[i,j] - sum_ldl)/d[j];

		end

		sum_dl = 0;

		for k=1:j

			sum_dl += d[k]*L[j+1,k]^2; 

		end


		d[j+1] = A[j+1,j+1] - sum_dl;

		 


	end
 

	return L, d;



end


#L_jul = F.L;
#D_jul = Diagonal(F.D);


# Example system Ax = b
A = [4.0 -1.0 0.0; -1.0 4.0 -1.0; 0.0 -1.0 3.0]
b = [15.0; 10.0; 10.0]
x0 = zeros(3)
k_max = 25

# Solve using the corrected Jacobi method
x = jacobi_method(A, b, x0, k_max)

println("Solution x: ", x)

println("Solution x*: ", A\b);
