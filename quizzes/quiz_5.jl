using Plots, LinearAlgebra
#part a

x = [1, 2, -0.5, -2 ,2];
y = [2, 3.35631, 0.76225, 0.126, 4.2563];

n = 5;
A =  [x'*x sum(x) ; sum(x) n+1];
b = [x'*y ; sum(y)];


Coeff = A \ b; 
display(Coeff);

#part b

log_y = log.(y);

A_exp = [x'*x sum(x); sum(x) n+1];
b_exp = [x'*log_y; sum(log_y)];


Coeff_exp = A_exp \ b_exp;
a_exp = Coeff_exp[1];
b_exp = exp(Coeff_exp[2]);


#part c
range = -3:0.01:3;
p(x) = @. Coeff[1]*x + Coeff[2];
q(x) = @. b_exp* exp(a_exp*x);
println("Exponential model coefficients: a  = ", a_exp, " b = " , b_exp);
scatter(x, y, label="Data Points", legend=:topleft);
plot!(range, p.(range), label="Linear Fit", title="Linear and Exponential Regression");
plot!(range, q.(range), label="Exponential Fit"); 
savefig("Plot.png");
#plot!(x, exp_fit.(x), label="Exponential Fit");

pts = [1,2, -0.5, -2 , 2];
#part d
res_linear = (y - p.(pts));
res_exp = (y - q.(pts));

max_res_lin = maximum((res_linear));
max_res_exp = maximum((res_exp));


cumulative_res_lin = sum(abs.(res_linear));
cumulative_res_exp = sum(abs.(res_exp));
println("Maximum absolute residual (Linear): ", max_res_lin);
println("Maximum absolute residual (Exponential): ", max_res_exp);

println("Cumulative residual (Linear): ", cumulative_res_lin);
println("Cumulative residual (Exponential): ", cumulative_res_exp);
