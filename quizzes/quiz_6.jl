using LinearAlgebra

A = [1 1 0;
     0 1 1;
     -1 1 -1]


U,S, V = svd(A)


A_inv = inv(A)
U_inv, S_inv, V_inc = svd(A_inv)


println("Singular values of A:", S);
println("Singular values of A^-1:", S_inv);


for i=1:3
	check = S[i]*S_inv[3-i+1];
	println(check);
end 

det_A = prod(S);
println("Determinant of A: ", det_A);
