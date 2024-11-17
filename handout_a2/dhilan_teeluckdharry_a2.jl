using LinearAlgebra



"""
Helper

"""
function gaussian_elimination(A, b)
    n = size(A,1);

    #Copying matrices to output variables
    Aout = copy(A);
    bout = copy(b);


    for k in 1:n-1
        for i in k+1:n

            #calculate the multiplier between pivot
            #and subsequent rows
            m_ik = Aout[i,k]/Aout[k,k];

            #we can set this to zero to 
            #reduce influence of cancellation
            Aout[i,k] = 0.0;

            #perform the row operation by substracting
            #row i by m_ik * row k 
            for j in k+1:n

               Aout[i,j] -= m_ik*Aout[k,j];

            end
            #update b_out as well
            bout[i] -= m_ik*bout[k];
        end
    end


    #return Aout and bout after GE
    return Aout, bout;

end


function backward_substitution(U, b)
    #defining solution vector x
    x = Vector{Float64}(undef, size(U,1));
    n = size(U,1);

    #for efficiency, x_n is calculated in O(1) here
    x[n] = b[n]/U[n,n];


    #back substitution algorithm
    for k in n-1:-1:1
        row_sum = 0;
        for i in n:-1:k+1
            #using subsequent values x_{k+1}, ... x_n
            #to calculate sum (u_i x_i ) for i > k
            row_sum += U[k,i]*x[i];

        end
        #calculating x_k
        x[k] = (b[k] - row_sum)/(U[k,k]);

    end

    #returning x
    return x;
end


"""
Compute the LDLᵀ decomposition (without pivoting) of a full-rank symmetric matrix A.
You may assume that A has a unique LU factorization. 

DO NOT COMPUTE THE LU DECOMPOSITION AS PART OF YOUR SOLUTION! Use the more efficient
procedure described in lecture.

Input:
    A: full-rank nxn symmetric real matrix with a unique LU factorization

Outputs:
    L: nxn real lower triangular matrix with identity on the diagonal
    d: vector of n real numbers such that A = L*diagm(d)*L'
"""
function ldl_decomposition(A)

	#initializing
	n = size(A,1);
        d = Vector{Float64}(undef,n);
        L = Matrix{Float64}(I,n,n);

        d[1] = A[1,1];
	
        for j = 1:n-1
		
		#computing each coefficent of L
                for i = j+1:n

                        sum_ldl = 0;

                        for k = 1:j-1
                                sum_ldl += L[i,k]*d[k]*L[j,k];
                        end

                        L[i,j] = (A[i,j] - sum_ldl)/d[j];

                end

                sum_dl = 0;

		#computing each coeffient of D

                for k=1:j

                        sum_dl += d[k]*L[j+1,k]^2;

                end


                d[j+1] = A[j+1,j+1] - sum_dl;




        end



    	return L, d
end

"""
Solve A*x = b for x given the LDL' = A decomposition for nxn symmetric matrix A.

Inputs:
    L: nxn lower triangular matrix
    d: n-element vector containing the diagonal entries of diagonal matrix D
    b: n-element vector (right-hand-side of the equation Ax=b)

Output:
    x: solution to A*x = b
"""
function ldl_solve(L, d, b)
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


	#init
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
	

	N = -(A-M);

        #G = M\N;
	
	G = zeros(Float64, n,n);

	for i =1:n
		for j=1:n
			if i!=j
				G[i,j] = N[i,j]/A[i,i];
			end
		end
	end

        cur_bound = norm(x_prev-x_true);

        push!(bounds, cur_bound);
	
	
	
        while iter < k_max && norm(A*x-b) >= res_tol

		#updating error

                push!(errs, norm(x_prev - x_true));
		
		#jacobi method for finding sol x
                for i = 1:n
                        sum_aij_xj = 0;

                        for j = 1:n
                                if j !== i
                                        sum_aij_xj += A[i,j]*x_prev[j];
                                end
                        end

                        x[i] = (b[i] - sum_aij_xj)/A[i,i];

                end

		#updating bound

                x_prev = copy(x);

                iter+=1;

                cur_bound = norm(G)*cur_bound;

                push!(bounds, cur_bound);

        end


        return x, errs, bounds;


end

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
function invert_lower(L)
	n = size(L,1);
	B = zeros(Float64,n,n);

	
	#loop over each row
	for i =1:n
		#diagonal element
		B[i,i] = 1/L[i,i];


		#loop over each column
		for j = 1:i-1
			sum = 0;
			for k = j:i-1
				sum += L[i,k]*B[k,j];
			end
			#compute off-diagonal element
			B[i,j] = -sum*B[i,i];

		end
		
		
	end

	return B;

end



function gauss_seidel(A, b, x0, x_true, k_max, res_tol)

	#init
	n = size(A,1);
        x_prev = copy(x0);
        iter  = 0;
        x = Vector{Float64}(undef,n);

        errs = Float64[];
        bounds = Float64[];

        N = zeros(Float64,n,n);
        M = Matrix{Float64}(I,n,n);

        #finding G matrix

        for i = 1:n
                for j=1:n
                        if j > i
                        N[i,j] = -A[i,j];
                        end
                end
        end

        M = N + A;

        #G = M\N;
	
	#G = zeros(Float64, n,n);
	
	M_inv = invert_lower(M);

	G = M_inv * N; 



        cur_bound = norm(x_prev-x_true);

        push!(bounds, cur_bound);



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

			#finding xi

                        x[i] = (b[i] - cur_elem_sum - cur_prev_sum)/(A[i,i]);




                end

		#updating error

                push!(errs, norm(x-x_true));

                x_prev = copy(x);
                iter+=1;


		#updating bound
                cur_bound = norm(G)*cur_bound;

                push!(bounds, cur_bound);


        end

        return x, errs, bounds;

end


"""
Computes the maximum magnitude eigenpair for the real symmetric matrix A 
with the power method. Ensures that the error tolerance is satisfied by 
using the Bauer-Fike theorem.

Inputs:
    A: real symmetric matrix
    tol: error tolerance; i.e., the maximum eigenvalue estimate 
         must be within tol of the true maximum eigenvalue

Outputs:
    λ: the estimate of the maximum magnitude eigenvalue
    v: the estimate of the eigenvector corresponding to λ
"""
function power_method_symmetric(A, tol)

	#init
        n = size(A,1);
        v = ones(n);

        lambda = 1;
        r = A*v - lambda*v;

        #for symmetric matrices, the bound evaluates to ||r||_2 /||x||_2
        bound = norm(r)/(norm(v));

	#applying power method until bound is below tol
        while bound >= tol

		#applying power method

                v_tilde = v/norm(v);
                v = A*v_tilde;

                lambda = v[1]/v_tilde[1];

		#updating bound

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
function extremal_eigenpairs(A, k, tol)

	#init
	n = size(A,1);
        lambda_vect = Vector{Float64}(undef,k);
        V = Matrix{Float64}(undef,n,k);

        A_mod = copy(A);

	#applying derived methodology

        for i=1:k
                lambda, eigvect = power_method_symmetric(A_mod,tol);
                lambda_vect[i] = lambda;

                V[:,i] = eigvect;



                A_mod -= lambda*eigvect*transpose(eigvect);

        end


        return lambda_vect, V;


end

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
	
	#helper function to compute jacobian at each iteration

        J = zeros(n,n);
	
	#using derived result

        for i = 1:n


                for j = 1:n

                        J[i,j] = (x[j] - p[j,i])/norm(x .- p[:,i]);

                end

        end

        return J;


end

function compute_fx(x,p,d,n)

	#helper function to compute fx at each iteration


        result = zeros(n);

	#using derived result

        for i = 1:n

                result[i] = norm(x .- p[:,i]) - d[i];

        end

        return result;

end



function newton(x0, P, d, tol, max_iters)

	#init

        n = size(P,1);
        x_trace = Vector{Vector{Float64}}();

        k = 0;

        x = copy(x0);

        fx_k = ones(n);


        #Jacobian

        while k <= max_iters && norm(fx_k) >= tol

                J = compute_jacobian(x,P,n);

                fx_k = compute_fx(x,P,d,n);


		#updating x
		U,b = gaussian_elimination(J,-fx_k);
		s = backward_substitution(U,b);


                #s = J \ -fx_k;

                x += s;

		#pushing to x_trace vector

                push!(x_trace, copy(x));

                k+=1;

        end


        return x_trace
end


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

	#helper function to compute grad f

        m = size(P)[2];
        n = size(P)[1];

        grad_f = zeros(n);

	#using derived result


        for i = 1:m
                diff = x .- P[:,i];
                dist = norm(diff);
                f_i = dist - d[i];

                grad_f .+= 2 * f_i * (diff/dist);


        end


        return grad_f;

end


function compute_hessian(x,P,d)

	#helper function to compute hessian

        m = size(P)[2];
        n = size(P)[1];

        grad_gradf = zeros(n,n);



	#using derived result
        for i = 1:m

                diff = x .- P[:,i];
                dist = norm(diff);
                f_i = dist - d[i];

                grad_gradf .+=  2*((diff * diff')/dist^2  + f_i * (I/dist  - (diff*diff')/(dist)^3) );

        end


        return grad_gradf;



end



function newton_optimizer(x0, P, d, tol, max_iters)

	#init

        m = size(P)[2];
        n = size(P)[1];

        x_trace = Vector{Vector{Float64}}();

        k = 0;


        x = copy(x0);

        grad_f = ones(n);

        while k <= max_iters && norm(grad_f) >= tol

                grad_f = compute_gradf(x,P,d);

                grad_gradf = compute_hessian(x,P,d);

		#updating x

		L,diag = ldl_decomposition(grad_gradf);
		
		s = ldl_solve(L,diag,-grad_f);

                s = grad_gradf \ -grad_f;

                x += s;

		#pushing to x_trace

                push!(x_trace, copy(x));


                k +=1;

        end

	#printing

	println("grad_f(x*) is:", compute_gradf(x, P,d));
	println("grad_gradf(x*) is:", compute_hessian(x,P,d));
	println("Eigenvalues of grad_gradf(x*) are:", eigvals(compute_hessian(x,P,d)));

        return x_trace

end
