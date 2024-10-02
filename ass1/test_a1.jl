using LinearAlgebra

# Replace with your solution:
include("template_a1.jl")


"""
Use this script to test your solutions (you will be graded by a similar script).
You may find it useful to modify the parameters or add your own tests.

DO NOT SUBMIT THIS FILE WITH YOUR SOLUTION!
"""

n_trials = 10
n = 100

## Test backward substitution
bs_score = 0
for _ in 1:n_trials
    A = rand(n, n)
    U = triu(A)
    b = rand(n)
    x_true = U \ b
    x_ans = backward_substitution(U, b)
    is_correct = all(abs.(x_true - x_ans) .<= sqrt(eps())*norm(x_true))
    global bs_score += is_correct
    if !is_correct
        display([x_true x_ans])
    end
end
println("Backward substitution score: $(bs_score)/$(n_trials)")

## Test forward substitution
fs_score = 0
for _ in 1:n_trials
    A = rand(n, n)
    L = tril(A)
    b = rand(n)
    x_true = L \ b
    x_ans = forward_substitution(L, b)
    is_correct = all(abs.(x_true - x_ans) .<= sqrt(eps())*norm(x_true))
    global fs_score += is_correct
    if !is_correct
        display([x_true x_ans])
    end
end
println("Forward substitution score: $(fs_score)/$(n_trials)")

## Test Gaussian elimination without pivoting (in conjunction with backward substitution)
ge_score = 0
for _ in 1:n_trials
    A = rand(n, n)
    b = rand(n)
    x_true = A \ b
    A_uef, b_uef = gaussian_elimination(A, b)
    x_ans = backward_substitution(A_uef, b_uef)
    is_correct = all(abs.(x_true - x_ans) .<= sqrt(eps())*norm(x_true))
    global ge_score += is_correct
    if !is_correct
        display([x_true x_ans])
    end
end
println("Gaussian elimination score: $(ge_score)/$(n_trials)")

## Test LU decomposition with partial pivoting
lu_score = 0 
for _ in 1:n_trials
    A = rand(n, n)
    b = rand(n)
    x_true = A \ b
    L, U, p = lu_partial_pivoting(A)
    x_ans = lu_solve(L, U, p, b)
    is_correct = all(abs.(x_true - x_ans) .<= sqrt(eps())*norm(x_true))
    global lu_score += is_correct
    if !is_correct
        display([x_true x_ans])
    end
end
println("LU decomposition score: $(lu_score)/$(n_trials)")
