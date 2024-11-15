using LinearAlgebra
using Plots



l0 = (x) -> @. x*(x-1)*(x-1.3)*(x-2)/(-1*-2*-2.3*-3);

l1 = (x) -> @. (x+1)*(x-1.3)*(x-1)*(x-2)/(1*-1.3*-1*-2);

l2 = (x) -> @. (x+1)*x*(x-1.3)*(x-2)/(2*1*-0.3*-1);

l3 = (x) -> @. (x+1)*x*(x-1)*(x-2)/(2.3*1.3*0.3*-0.7);

l4 = (x) -> @. (x+1)*x*(x-1)*(x-1.3)/(3*2*1*0.7);



p = (x) -> @. -1*l0(x) + 2*l1(x) + 1 * l2(x) + 0.967*l3(x) + 3*l4(x);


x = -1:0.1:2;

plot(x,p(x));
savefig("lagrange.png");

scatter!([-1,0,1,1.3,2], [-1,2,1,0.967,3], legend=false);


f = (x) -> @. 7/6*x^3 - 2*x^2 - 1/6*x + 2;

max_abs_err = 0;

g_vals = p(x);
f_vals = f(x);

max_abs_err = maximum(abs.(f_vals .- g_vals));

println(max_abs_err);





