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

function extremal_eigenpairs(A,k,tol)
	n = size(A,1);
	lambda_vect = Vector{Float64} (undef,n);

	


end



