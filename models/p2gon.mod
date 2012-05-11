# Maximum area for unit-diameter polygon (Prieto's second model).

# This derives from a GAMS model by Francisco J. Prieto -- the second
# one in E-mail that he sent to Margaret Wright.  Prieto described
# this as "a very particular solution for the problem, taking advantage
# of [a] certain property that I [Prieto] am convinced (although I cannot
# prove it) that the solution must have.  The formulation
# is much more efficient, and you can go up to sizes of around 70-80.

# "Symmetry is assumed (so the size is reduced to a half), and also as
# many constraints as possible are set up as equalities.

# "The results obtained seem to coincide (within some error margin) with
# the ones obtained from the previous program."


param N 'number of sides' integer > 1;
	check: N == 2*floor(N/2);	# Make sure N is even.

param Lsup := N/2 - 1;

set I := 1..Lsup;

var x{I} >= 0, <= 1, := .5;		#first (Cartesian) coordinate

var y{i in I} >= 0, <= 1, := i/Lsup;	#second (Cartesian) coordinate

# distance constraints (tight)
s.t. cd{i in ceil(Lsup/2) .. Lsup, j in max(1, Lsup-i) .. min(Lsup+1-i, i)}:
	(x[i] + x[j])**2 + (y[i] - y[j])**2 == 1;

# extra distance constraint
s.t. extr: x[Lsup]^2 + y[Lsup]^2 == 1;

# second coordinate constraint (ordered values)
s.t. ory{i in 2..Lsup}: y[i] >= y[i-1];

maximize area: sum{i in 2..Lsup-1} y[i]*(x[i-1]-x[i+1])  +
		x[Lsup] + y[Lsup]*x[Lsup-1] - x[2]*y[1];

data;
param N :=

# This assumes that N will be given separately.  For example,
# you might invoke
#    ampl -omp20 p2gon.mod -
# and then type
#    20; end;
# to generate the 20-sided version of this problem (in files p20.*).
