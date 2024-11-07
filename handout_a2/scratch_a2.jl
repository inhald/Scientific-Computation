using LinearAlgebra


"""
Solve Ax=b iteratively using the Gauss-Seidel method with initial iteration x0.

Inputs:
    A: nxn real matrix
    b: n-element real vector
    x0: initial iteration
    x_true: ground-truth value used to compute output errs
    k_max: maximum number of iterations
    res_tol: residual error tolerance criterion; terminate once norm(A*x - b) <= res_tol
Outputs:
    x: n-element real vector containing the solution returned by the Gauss-Seidel method
    errs: vector of k+1 Euclidean norm error values, where k is the number of
          iterations performed. The i-th element errs[i] contains the absolute error
          in x after i-1 iterations of the Gauss-Seidel method have been performed (e.g.,
          errs[1] contains norm(x0 - x_true)).
    bounds: upper bound on errs computed using the convergence theory discussed in lecture.
"""
function gauss_seidel(A, b, x0, x_true, k_max, res_tol)

	n = size(A,1);
	x_prev = copy(x0);
	iter  = 0;
	x = Vector{Float64}(undef,n);

	err = Float64[];
	bounds = Float64[];




	while iter <= k_max && norm(A*x - b) >= res_tol


		for i=1:n
			#cur sum
			cur_elem_sum = 0;

			for j =1:i-1
				cur_elem_sum += A[i,j] * x[j];

			end

			#prev sum 
			cur_prev_sum = 0;


			for j=i+1:n
				cur_prev_sum += A[i,j]*x_prev[j];			

			end

			x[i] = (b[i] - cur_elem_sum - cur_prev_sum)/(A[i,i]);


			

		end

		push!(err, norm(x-x_prev));

		x_prev = copy(x);
		iter+=1; 




	end



    	return x
end


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
function jacobi_method(A, b, x0, x_true, k_max, res_tol)

	n = size(A,1);
	x = Vector{Float64}(undef,n);
	errs = Float64[];
	x_prev = copy(x0);

	iter = 0;
	

	bounds = Float64[];
	M = Matrix{Float64}(I,n,n);

	#finding G matrix
	
	for i = 1:n
		M[i,i] = A[i,i];
	end

	N = M - A;  

	G = M\N;

	cur_bound = x_prev-x_true;
	
	push!(bounds, norm(cur_bound));

	
	while iter <= k_max && norm(A*x-b) >= res_tol

		push!(errs, norm(x_prev - x_true));

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

		iter+=1; 
		
		cur_bound = G*cur_bound;

		push!(bounds, norm(cur_bound));

	end


    	return x, errs, bounds; 
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


end


#L_jul = F.L;
#D_jul = Diagonal(F.D);


# Example system Ax = b
A = [4.0 -1.0 0.0; -1.0 4.0 -1.0; 0.0 -1.0 3.0]
b = [15.0; 10.0; 10.0]
x0 = zeros(3)
k_max = 25

# Solve using the corrected Jacobi method
x, err, bounds = jacobi_method(A, b, x0, A\b,  k_max, 1e-6);

x = gauss_seidel(A,b,x0, A\b, k_max, 1e-6);

iter = 0; 

iter+=1; 

#println(iter);


println("Solution x: ", x)

println("Solution x*: ", A\b);
