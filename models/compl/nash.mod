#==>nash.gms

# A non-cooperative game example: a Nash equilibrium is sought.

# References:

# F.H. Murphy, H.D. Sherali, and A.L. Soyster, "A Mathematical
# 	Programming Approach for Determining Oligopolistic Market
# 	Equilibrium", Mathematical Programming 24 (1986), pp. 92-106.

# P.T. Harker, "Accelerating the convergence . . .", Mathematical
# 	Programming 41 (1988), pp. 29-59.

set Rn := 1 .. 10;

param gamma := 1.2;

param c {i in Rn};
param beta {i in Rn};
param L {i in Rn} := 10;

var q {i in Rn} >= 0;		# production vector
var Q = sum {i in Rn} q[i];
var divQ = (5000/Q)**(1/gamma);

s.t. feas {i in Rn}:
	q[i] >= 0 complements
	0 <= c[i] + (L[i] * q[i])**(1/beta[i]) - divQ
	     - q[i] * (-1/gamma) * divQ / Q;
