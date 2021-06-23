# gridneta.mod	OBR2-MY-V60-V96
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs);
# modified to permit variable problem sizes and bound choices by David M. Gay.

# Ref.: Ph.Toint and D.Tuyttens, On Large Scale Nonlinear Network Optimization,
# Mathematical Programming B, vol.48(1), pp. 125-159, 1990.

# Params L, a, c, ibounds, and r give what the ref. calls
# problem P(L,a,c,ibounds,r).
# The default settings give P(2,1,2,0,.1).

# Number of variables: 60 (with default settings)
# Number of constraints: 96 (with default settings)
# Objective partially separable
# Objective convex (?)
# Linear constraints


param L default 2;
param M := 2*L + 1;
param nv := 2*M*(M+1);	# number of variables (for complicated bound spec.)
set I  := 0 .. M;
set I1 := 0 .. M-1;
node N{i in I, j in I}:
	net_out == if i == 0 && j == 0 then +10 else
		   if i == M && j == M then -10 else 0;

set A := {i in I, j in I, k in 0..1: (i < M || k == 1) && (j < M || k == 0)};

arc x{(i,j,k) in A} >= 0 from N[i,j] to N[i+1-k,j+k];

param J{(i,j,k) in A} :=
	if k == 0 then (if j > 0 && j < M then floor((j+1)/2)
			else if i > 0 && i < M then (i mod 2)*(i+1)/2)
		  else (if i > 0 && i < M then floor((i+1)/2)
			else if j > 0 && j < M then (j mod 2)*(j+1)/2);
param Jprev{(i,j,k) in A} :=
		if k == 1 then (if i < M then J[i,j,0]
				else J[M-1,j,1])
			  else if i > 0 then (if j == M then J[i-1,M,0]
						else J[i-1,j,1])
			  else if j > 0 then J[M,j-1,1];

# Number facilitates the computatation of lower and upper bounds
# specified in the reference above.

param Number{(i,j,k) in A} := if max(i,j,k) == 0 then 1
			  else 1 + (if k == 1 then (if i < M then Number[i,j,0]
						else Number[M-1,j,1])
			  	else if i > 0 then
					(if j == M then Number[i-1,M,0]
						else Number[i-1,j,1])
			  	else Number[M,j-1,1]);

param c default 2;
param a default 1;
param r default .1;
param alpha {(i,j,k) in A}
		:= if J[i,j,k] >= 1 then 10^((Jprev[i,j,k]/(L-1))*log10(c))
				    else 1;

param sign{(i,j,k) in A} := (-1) ^ if k == 0 then (if j < M then j+1 else i)
					     else (if i < M then j else j+1);

minimize Cost:
          (1/100) * sum {(i,j,k) in A} alpha[i,j,k]*x[i,j,k]^2
	+ (1/100) * a
	  * ( sum {(i,j,k) in A: j < M || i < M-1}
		 sqrt(1 + x[i,j,k]^2
			+ (x[i,j,k] - if k == 0
				then (if j == M then x[i+1,M,0]
						else x[i,j,1])
				else (if i == M then x[0,j+1,0]
					else if i == M-1 then x[M,j,1]
					else x[i+1,j,0])
				)^2)
	    + (1/1200)*(10 + sum {(i,j,k) in A} sign[i,j,k]*x[i,j,k] )^4);

param ibounds integer in -1..1 default 0;	# bounds case

set Border :=	1 .. 2*M-1 by 2 union
		2*M+1 .. M*(2*M+1) by 2*M+1 union
		2 .. 2 + (M-1)*(2*M+1) by 2*M+1 union
		M*(2*M+1)+1 .. M*(2*M+2);

s.t. bounds{(i,j,k) in A, n in {Number[i,j,k]}}:
	if ibounds == -1 then
		if n mod 3 == 0 && n not in Border
		then -Infinity
		else if n in {nv - 3*L - 8, nv - L - 3, nv - 1, nv}
		then 4
		else if n in Border
		then r
		else -Infinity
	else if ibounds == 0 then
		-Infinity
	else if n mod 3 == 0 then r
	else -Infinity
	<= x[i,j,k]
	<=
	if ibounds == -1 then
		if n mod 3 == 0 && n not in Border
		then 0
		else if n in {nv - 3*L - 8, nv - L - 3, nv - 1, nv}
		then 5
		else Infinity
	else Infinity;
