import LinearAlgebra as linalg



import LinearAlgebra as linalg

function power_method(A, x, k_max)
    # Normalize initial vector x
    x /= linalg.norm(x)

    for k in 1:k_max
        # Multiply A by the current estimate of x
        y = A * x

        # Update x by normalizing y
        x = y / linalg.norm(y)

        # Approximate the dominant eigenvalue
        λ = linalg.dot(x, A * x)

        println("Iteration $k: Approximate eigenvalue = $λ, Eigenvector = $x")
    end
end

# Define the matrix A and initial vector x
A = [3 1; 1 3]
k_max = 10
x = [-0.5, 0.5]

# Run the Power Method
power_method(A, x, k_max)
