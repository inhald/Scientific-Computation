using Random
include("dhilan_teeluckdharry_a2.jl")

Random.seed!(987654321)

n = 200
A = rand(n, n)
α = 2.0
A = A+A' + α*n*I
x_true = rand(n)
b = A*x_true
x0 = rand(n)
res_tol = 0.0
k_max = 50

# Run the iterative solvers
x_jac, errs_jac, bounds_jac = jacobi_method(A, b, x0, x_true, k_max, res_tol)
jac_iters = 1:length(errs_jac)
x_gs, errs_gs, bounds_gs = gauss_seidel(A, b, x0, x_true, k_max, res_tol)
gs_iters = 1:length(errs_gs)
println(eps(x_jac[1]));

L, d = ldl_decomposition(A)
x_ldl = ldl_solve(L, d, b)
res_ldl = norm(A*x_ldl - b)

# Plot the results
using Plots

plot(errs_jac, yaxis=:log, label="Jacobi", color=:red)
plot!(bounds_jac, yaxis=:log, label="Jacobi Bound", color=:red, style=:dash)
plot!(errs_gs, yaxis=:log, label="Gauss-Seidel", color=:green)
plot!(bounds_gs, yaxis=:log, label="Gauss-Seidel Bound", color=:green, style=:dash)
hline!([res_ldl], yaxis=:log, label="LDL' Reference", color=:blue)
title!("Convergence of Iterative Methods")
xaxis!("Iteration")
yaxis!("Absolute Error")


savefig("Results.png");
