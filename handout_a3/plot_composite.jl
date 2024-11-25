using Plots

# Include your completed code
include("firstname_lastname_a3.jl")

# Define the function and bounds for our integral
f(x) = exp(-x/2)*sin(x)
a = 0.0
b = 4*π

# Compute the analytical solution
f_int(x) = -2*exp(-x/2)*(sin(x) + 2*cos(x))/5
true_int = f_int(b) - f_int(a)

# Compute approximate integrals for a variety of values r
r = 2 .^(2:20)
h = (b-a)./r
int_trap = [composite_trapezoidal_rule(f, a, b, r[i]) for i in eachindex(r)]
int_mid = [composite_midpoint_rule(f, a, b, r[i]) for i in eachindex(r)]
int_simpson = [composite_simpsons_rule(f, a, b, r[i]) for i in eachindex(r)]

# Plot the results on log axes
plot(h, abs.(int_trap .- true_int), label="Trapezoid", marker=:o, 
     yaxis=:log10, xaxis=:log10, legend=:topleft)
plot!(h, abs.(int_mid .- true_int), label="Midpoint", marker=:d)
plot!(h, abs.(int_simpson .- true_int), label="Simpson", marker=:s)
title!("Error of Composite Quadrature Rules")
xaxis!("h")
yaxis!("error")