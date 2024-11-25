using Plots

# Include your completed code
include("firstname_lastname_a3.jl")

# Problem parameters
ω = 2.0
tol = 0.2
a = -π
b = π
N = 1000 # For plotting
x = range(a, stop=b, length=N)

# Use your solution to subdivide (evenly spaced nodes)
n = subdivide(a, b, ω, tol)
x_data = range(a, stop=b, length=n)
y_data = cos.(ω*x_data)
c_even = newton_int(x_data, y_data)
p_even = horner(c_even, x_data, x)

# Use your solution to chebyshev_nodes
x_cheb = chebyshev_nodes(a, b, ω, tol)
y_cheb = cos.(ω*x_cheb)
c_cheb = newton_int(x_cheb, y_cheb)
p_cheb = horner(c_cheb, x_cheb, x)

plot(x, abs.(p_even .- cos.(ω*x)), label="Evenly-spaced Interpolation Error", color=:red)
scatter!(x_data, zeros(length(x_data)), label="Evenly-spaced Nodes", markercolor=:red)
plot!(x, abs.(p_cheb .- cos.(ω*x)), label="Chebyshev Interpolation Error", color=:green)
scatter!(x_cheb, zeros(length(x_cheb)), label="Chebyshev Nodes", markercolor=:green)
hline!([tol], label="Tolerance", color=:black, linestyle=:dash)