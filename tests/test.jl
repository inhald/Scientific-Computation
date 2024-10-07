# Quiz1 Attempt
using LinearAlgebra
import Random

function lu_partial_pivoting(A)
    n = size(A,1);
    L = Matrix{Float64}(I,n,n);
    U = Matrix{Float64}(undef,n,n);

    p = Vector{Float64}(undef, n);
    
    for i=1:n
        p[i] = i;
    end

    for k = 1:n-1
        #find pivot in each column
        pivot = k;
        col_max = abs(A[k,k]);
                
        for i = k+1:n

            if abs(A[i,k]) >= col_max
                col_max = abs(A[i,k]);
                pivot = i;
            end

        end

        p[k], p[pivot] = p[pivot], p[k];

        # println(string(col_max));


        #exchanging rows
        if pivot != k

            for i = 1:n
                A[pivot, i], A[k,i] = A[k,i], A[pivot, i];
                L[pivot, i], L[k,i] = L[k,i], L[pivot, i];
            end
        end

        # println(string(A));


        #applying normal GE


        for i = k+1:n 

            m_ik = A[i,k]/A[k,k];
            L[i,k] = m_ik;
            A[i,k] = 0.0;
            
            for j=k+1:n
                A[i,j] -=  m_ik*A[k,j];
            end

        end
        

    end

    U = A;
    # println(string(A));

    return L, U, p;
end


A = [1.0 2.0 3.0; 3.0 4.0 5.0; 13.0 12.0 -19.0];


L,U,p = lu_partial_pivoting(A);

println(string(L));
println(string(U));
println(string(p));
