# Maximum area for unit-diameter polygon of N sides.
# The following model started as a GAMS model by Francisco J. Prieto.

param N integer > 0;
set I := 1..N;

param pi := 4*atan(1.);

var rho{i in I} <= 1, >= 0	# polar radius (distance to fixed vertex)
		:=  4*i*(N + 1 - i)/(N+1)**2;

var the{i in I} >= 0		# polar angle (measured from fixed direction)
		:= pi*i/N;

s.t. cd{i in I, j in i+1 .. N}:
	rho[i]**2 + rho[j]**2 - 2*rho[i]*rho[j]*cos(the[j]-the[i]) <= 1;

s.t. ac{i in 2..N}:
	the[i] >= the[i-1];

s.t. fix_theta: the[N] = pi;
s.t. fix_rho:   rho[N] = 0;

maximize area: .5*sum{i in 2..N} rho[i]*rho[i-1]*sin(the[i]-the[i-1]);

data;
param N :=

# This assumes that N will be given separately.  For example,
# you might invoke
#    ampl -omp6 pgon.mod -
# and then type
#    6; end;
# to generate the 6-sided version of this problem (in files p6.*).
