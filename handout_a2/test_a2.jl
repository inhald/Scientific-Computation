using Random
using LinearAlgebra
Random.seed!(99991099910101010)

# Replace with your solution:
include("template_a2.jl")

"""
Use this script to test your solutions (you will be graded by a similar script).
You may find it useful to modify the parameters or add your own tests.

DO NOT SUBMIT THIS FILE WITH YOUR SOLUTION!
"""

n_trials = 10
n = 100

"""
A test implementation of the LU decomposition without pivoting. 
Your LDL decomposition should be significantly faster than this! 
"""
function lu_decomposition(A)
    n = size(A, 1)
    A_out = copy(A)
    for k in 1:n-1
        inds = k+1:n
        # Perform a step of Gaussian elimination; store non-diagonal part of L and all of U in A!
        A_out[inds, k] = A_out[inds, k]/A_out[k, k]  # This step needs to go first (computing terms in L)
        A_out[inds, inds] = A_out[inds, inds] - A_out[inds, k]*A_out[k, inds]'  # This updates U and the rest of L
    end 
    return tril(A_out, -1) + I, triu(A_out)
end

## Test ldl_decomposition() accuracy
ldl_acc_score = 0
for _ in 1:n_trials
    A = 2.0*rand(n, n) .- 1.0
    A = A + A' + 10.0*I  # Encourage diagonal dominance
    L, d = ldl_decomposition(A)
    L_test, U_test = lu_decomposition(A)
    L_is_correct = all(abs.(L - L_test) .<= sqrt(eps()))
    d_is_correct = all(abs.(diagm(d)*L' - U_test) .<= sqrt(eps()))
    global ldl_acc_score += L_is_correct && d_is_correct
end
println("LDL accuracy score: $(ldl_acc_score)/$(n_trials)")

## Test ldl_decomposition() speed
ldl_speed_score = 0
n_time_trials = 50
for i in 1:n_time_trials+1
    A = rand(n, n)
    A = A + A' + 10.0*I
    b = rand(n)
    t_ldl = @elapsed ldl_decomposition(A)
    t_lu = @elapsed lu_decomposition(A)
    if i > 1 # Avoid compilation time
        global ldl_speed_score += t_ldl < 0.6 * t_lu
    end
end
ldl_speed_marks = Int(round(Float64(ldl_speed_score)/10))
println("LDL speed score: $(ldl_speed_marks)/5")
println("Problem 1 score: $(ldl_acc_score + ldl_speed_marks)/$(15) marks")

## Test ldl_solve() and ldl_decomposition() accuracy together
lu_score = 0 
for _ in 1:n_trials
    A = rand(n, n)
    A = A + A' + 10.0*I
    b = rand(n)
    x_true = A \ b
    L, d = ldl_decomposition(A)
    x_ans = ldl_solve(L, d, b)
    is_correct = all(abs.(x_true - x_ans) .<= sqrt(eps())*norm(x_true))
    global lu_score += is_correct
    if !is_correct
        display([x_true x_ans])
    end
end
println("Problem 2 score: $(lu_score)/$(n_trials) marks")

## Test power_method_symmetric accuracy
power_method_score = 0
tol = 1e-9
n_trials=50
for _ in 1:n_trials
    A = 2.0*rand(n, n) .- 1.0
    A = A + A'
    λ, v = power_method_symmetric(A, tol)
    eigpairs = eigen(A)
    _, max_ind = findmax(abs.(eigpairs.values))
    val_correct = abs(eigpairs.values[max_ind] - λ) <= tol 
    vec_correct = norm(eigpairs.vectors[:, max_ind] - v) <= tol*sqrt(n) ||
                  norm(eigpairs.vectors[:, max_ind] + v) <= tol*sqrt(n)
    global power_method_score += val_correct && vec_correct
end
power_method_score = Int(round(power_method_score/10.0))
println("Problem 4 score: $(power_method_score)/5 marks")

## Test extremal_eigenpairs accuracy
extremal_eigs_score = 0
for _ in 1:n_trials
    A = 2.0*rand(n, n) .- 1.0
    A = A + A'
    k = rand(2:5)
    λ, _ = extremal_eigenpairs(A, k, tol)
    eigpairs = eigen(A)
    max_inds = sortperm(abs.(eigpairs.values), rev=true)
    λ_true = eigpairs.values[max_inds[1:k]]
    vals_correct = all(abs.(λ_true - λ) .<= tol)
    global extremal_eigs_score += vals_correct
end
extremal_eigs_score = Int(round(extremal_eigs_score/10.0))
println("Problem 6 score: $(extremal_eigs_score)/5 marks")


# Test out the Newton solver 
score_newton = 0
newton_tol = 1e-6
for _ in 1:n_trials
    x_gt = rand(3)*10 .- 5.0 
    x0 = x_gt + 0.5*rand(3)
    P = rand(3, 3)*10.0 .- 5.0
    d = [norm(x_gt - P[:, i]) for i in 1:size(P, 2)]
    x_trace = newton(x0, P, d, newton_tol, 100)
    global score_newton += (norm(x_trace[end] - x_gt) <= newton_tol)
end
score_newton = Int(round(score_newton/5.0))
println("Problem 8 score: $(score_newton)/10")
