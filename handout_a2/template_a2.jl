using LinearAlgebra

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

        push!(bounds, norm(x_prev-x_true));


        while iter <= k_max && norm(A*x-b) <= res_tol

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

                push!(bounds, G*bounds[iter]);

        end



    	return x, errs, bounds
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
function gauss_seidel(A, b, x0, x_true, k_max, res_tol)
    return x, errs, bounds
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
    return λ, v
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
    return λ, V
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
function newton(x0, P, d, tol, max_iters)
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
function newton_optimizer(x0, P, d, tol, max_iters)
    return x_trace
end
