using LinearAlgebra
A = [1 0 1/3 0; 0 1 3 -1; 0 0 8 3; 0 0 0 -13/4]

Fact = LinearAlgebra.lu(A, LinearAlgebra.NoPivot(););

println(string(Fact.U));

