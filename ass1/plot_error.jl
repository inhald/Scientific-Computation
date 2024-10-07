using LinearAlgebra
using Statistics
using Plots

"""
Use this script to generate the plots required for Problem 8 - be sure to modify 
the savefig function call at the end of this file to save your figures.

DO NOT SUBMIT THIS FILE WITH YOUR SOLUTION!
"""

# Replace with your solution:
include("template_a1.jl")

n = 200
n_trials = 200

res_ge = zeros(n_trials)
err_ge = zeros(n_trials)
res_lu = zeros(n_trials)
err_lu = zeros(n_trials)
condition_numbers = zeros(n_trials)

for trial in 1:n_trials
    A = rand(n, n)
    x = rand(n)
    b = A*x
    A_ge, b_ge = gaussian_elimination(A, b)
    x_ge = backward_substitution(A_ge, b_ge)
    L, U, p = lu_partial_pivoting(A)
    x_lu = lu_solve(L, U, p, b)
    err_ge[trial] = norm(x_ge - x)
    err_lu[trial] = norm(x_lu - x)
    res_ge[trial] = norm(A*x_ge - b)
    res_lu[trial] = norm(A*x_lu - b)
    condition_numbers[trial] = cond(A)
end

p1 = scatter(condition_numbers, err_ge, color=:red, label="Gaussian Elimination", xaxis=:log, yaxis=:log)
scatter!(condition_numbers, err_lu, color=:green, label="Partial Pivoting", xaxis=:log, yaxis=:log)
xaxis!("Condition Number")
yaxis!("Euclidean Norm of Solution Error")
title!("Error vs. Condition Number")
display(p1)

p2 = scatter(condition_numbers, res_ge, color=:red, label="Gaussian Elimination", xaxis=:log, yaxis=:log)
scatter!(condition_numbers, res_lu, color=:green, label="Partial Pivoting", xaxis=:log, yaxis=:log)
xaxis!("Condition Number")
yaxis!("Euclidean Norm of Residual")
title!("Residual vs. Condition Number")
display(p2)


n_trials_α = 100
α_diag = 10.0 .^ (-8:1:8)
n_α = length(α_diag)
res_ge_α = zeros(n_α, n_trials_α)
res_lu_α = zeros(n_α, n_trials_α)
err_ge_α = zeros(n_α, n_trials_α)
err_lu_α = zeros(n_α, n_trials_α)

for (i, α) in enumerate(α_diag)
    for trial in 1:n_trials_α
        A = rand(n, n) + n*α*I
        x = rand(n)
        b = A*x
        A_ge, b_ge = gaussian_elimination(A, b)
        x_ge = backward_substitution(A_ge, b_ge)
        L, U, p = lu_partial_pivoting(A)
        x_lu = lu_solve(L, U, p, b)

        err_ge_α[i, trial] = norm(x_ge - x)
        err_lu_α[i, trial] = norm(x_lu - x)        
        res_ge_α[i, trial] = norm(A*x_ge - b)
        res_lu_α[i, trial] = norm(A*x_lu - b)        
    end
end

mean_res_ge = mean(res_ge_α, dims=2)
mean_res_lu = mean(res_lu_α, dims=2)
σ_res_ge = std(res_ge_α, dims=2)
σ_res_lu = std(res_lu_α, dims=2)
p_mean_res = plot(α_diag, mean_res_ge, label="Gaussian Elimination", color=:red, xaxis=:log, yaxis=:log)
plot!(α_diag, mean_res_lu, label="Partial Pivoting", color=:green, xaxis=:log, yaxis=:log, linestyle=:dash)
xaxis!("Diagonal Weight α")
yaxis!("Mean Euclidean Norm of Residual")
title!("Residual vs. Diagonal Scaling ($(n_trials_α) Trials per α)")
display(p_mean_res)

mean_err_ge = mean(err_ge_α, dims=2)
mean_err_lu = mean(err_lu_α, dims=2)
p_mean_err = plot(α_diag, mean_err_ge, label="Gaussian Elimination", color=:red, xaxis=:log, yaxis=:log)
plot!(α_diag, mean_err_lu, label="Partial Pivoting", color=:green, xaxis=:log, yaxis=:log, linestyle=:dash)
xaxis!("Diagonal Weight α")
yaxis!("Mean Euclidean Norm of Error")
title!("Error vs. Diagonal Scaling ($(n_trials_α) Trials per α)")
display(p_mean_err)

# Save the relevant files to add them to your PDF solution, e.g.:
savefig(p1, "p1.png");
savefig(p2, "p2.png");
savefig(p_mean_err, "p_mean_err.png");
savefig(p_mean_res, "p_mean_res.png" );