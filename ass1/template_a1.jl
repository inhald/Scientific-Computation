using LinearAlgebra

"""
Implement the function signatures defined here in your solution file <FIRSTNAME>_LASTNAME_a1.jl.
DO NOT MODIFY THE SIGNATURES.

Use test_a1.jl to test your solutions.
"""


"""
Perform backward substitution to solve the system U*x = b for x.
Inputs:
    U: an nxn upper diagonal square matrix (assume the diagonal is nonzero)
    b: an nx1 vector
Output: 
    x: the solution to U*x = b
"""
function backward_substitution(U, b)
    x = Vector{Float64}(undef, size(U,1));
    n = size(U,1);

    x[n] = b[n]/U[n,n];


    for k in n-1:-1:1
        my_sum = 0;
        for i in n:-1:k
            my_sum += U[k,i]*x[i];

        end

        x[k] = (b[k] - my_sum)/(U[k,k]);

    end

    return x;
end

"""
Perform forward substitution to solve the system L*x = b for x.
Inputs:
    L: an nxn lower diagonal square matrix (assume the diagonal is nonzero)
    b: an nx1 vector
Output: 
    x: the solution to L*x = b
"""
function forward_substitution(L, b)
    x = Vector{Float64}(undef, size(L,1));
    n = size(L,1);

    x[1] = b[1]/L[1,1];

    for k in 2:n
        elem_sum = 0;
        for i = 1:k-1
            elem_sum += L[k,i]*x[i];

        end

        x[k] = (b[k] - elem_sum)/(L[k,k]);

    end

    return x; 
end

"""
Gaussian elimination without pivoting. 
You may assume A_in has full rank and never has a zero pivot. 

Input:
    A: a full rank nxn matrix with no zero pivots.
    b: a nx1 vector

Output:
    A_out: the input matrix A in row-echelon form
    b_out: the input vector b after having undergone the transformations putting A into row-echelon form
"""
function gaussian_elimination(A, b)
    n = size(A,1);

    Aout = copy(A);
    bout = copy(b);

    for k in 1:n-1
        for i in k+1:n

            if Aout[i,k] < eps(Aout[i,k])
            end

            m_ik = Aout[i,k]/Aout[k,k];
            
            Aout[i,k] = 0.0; 
            
            for j in k+1:n

               Aout[i,j] -= m_ik*Aout[k,j];

            end
            
            bout[i] -= m_ik*bout[k];
        end
    end


    return Aout, bout;

end

"""
LU decomposition with partial pivoting. 

Input: 
    A_in: an nxn full-rank matrix

Returns:
    L: nxn lower triangular matrix
    U: nxn upper triangular matrix
    p: permuted vector of indices 1:n representing the pivots 
       (i.e., your solution should (approximately) satisfy A[p, :] = L*U) 
"""
function lu_partial_pivoting(A)
    n = size(A,1);
    L = Matrix{Float64}(I,n,n);
    p = Vector{Int}(undef, n);
    U = copy(A);
    
    for i=1:n
        p[i] = i;
    end

    for k = 1:n-1
        #find pivot in each column
        pivot = k;
        col_max = abs(U[k,k]);
                
        for i = k+1:n

            if abs(U[i,k]) >= col_max
                col_max = abs(U[i,k]);
                pivot = i;
            end

        end


        #exchanging rows
        if pivot != k

            # for i = 1:n
            #     U[pivot, i], U[k,i] = U[k,i], U[pivot, i];
            # end
            U[[pivot,k], 1:n] = U[[k, pivot], 1:n];

            p[k], p[pivot] = p[pivot], p[k];


            if k >= 2
                # for i = 1:k-1
                    # L[k, i] , L[pivot, i] = L[pivot, i] , L[k, i];
                # end

                # L[[k,pivot]: 1:k-1] = L[[pivot,k], 1:k-1];
                L[[k, pivot], 1:k-1] = L[[pivot, k], 1:k-1];

            end

        end

        #applying normal GE


        for i = k+1:n 

            m_ik = U[i,k]/U[k,k];
            L[i,k] = m_ik;
            
            for j=k+1:n
                U[i,j] -=  m_ik*U[k,j];
            end
            U[i,k] = 0.0;

        end
        

    end

    return L, U, p; 
 
end


"""
Solve A*x = b for x given the LU decomposition (approximately) satisfying A[p, :] = L*U.

Inputs:
    L: lower triangular matrix
    U: upper triangular matrix
    p: vector of permuted indices representing the pivots

Output:
    x: solution to A*x = b (recall that A[p, :] = L*U)
"""
function lu_solve(L, U, p, b)
    
    n = size(L,1);
    b_prime = Vector{Float64}(undef, n);

    for i = 1:n
        b_prime[i] = b[p[i]];
    end

    z = forward_substitution(L,b_prime);

    x = backward_substitution(U,z);

    return x;
end
