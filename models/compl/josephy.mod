# josephy.mod

# Complementarity problem from Kojima, via Josephy.

set Rn;

param A {Rn,Rn,Rn};
param B {Rn,Rn};
param c {Rn};

var x {j in Rn};

s.t. f {i in Rn}:
	0 <= x[i]
     complements
	sum {j in Rn} (B[i,j]*x[j] + sum {k in Rn} A[i,j,k]*x[j]*x[k])
	>= -c[i];
