#  obstacle.mod

# Calculate the position of a membrane pushed up through a (rectangular)
# hole in a rigid plate; in addition, there are rigid obstacle(s) inside
# the hole (perhaps not at the same level as the plate).  These
# obstacles can constrain the membrane from above or from below.  The
# correct position of the membrane (the position of minimum energy) is
# determined by the minimization of a quadratic function of the membrane
# position, subject to the constraints imposed by the hole and the
# obstacle.  The MCP below arises from the optimality conditions for
# this QP.

# Reference: Ciarlet, Philippe G., "The Finite Element Method for
# Elliptic Problems", North-Holland, 1978.

param M integer >= 1, default 50; # of interior grid pts in Y direction
param N integer >= 1, default 50; # of interior grid pts in X direction

set Y := 0 .. M+1 ;
set X := 0 .. N+1 ;

param xlo := 0;
param xhi > xlo, := 1;
param ylo := 0;
param yhi > ylo, := 1;
param dy := (yhi - ylo) / (M + 1) ;
param dx := (xhi - xlo) / (N + 1) ;
param c := 1;		/* force constant */

param ub {i in Y, j in X} := if 1 <= i <= M and 1 <= j <= N then
			(sin(9.2*(xlo+dx*i))*sin(9.3*(ylo+j*dy)))^2 + 0.2;
param lb {i in Y, j in X} := if 1 <= i <= M and 1 <= j <= N then
			(sin(9.2*(xlo+dx*i))*sin(9.3*(ylo+j*dy)))^3;

# height of membrane
var v {i in Y, j in X} >= lb[i,j]  <= ub[i,j]  := max(0,lb[i,j]);

s.t. dv {i in 1..M, j in 1..N}:
	lb[i,j] <= v[i,j] <= ub[i,j] complements
		  (dy/dx) * (2*v[i,j] - v[i+1,j] - v[i-1,j])
		+ (dx/dy) * (2*v[i,j] - v[i,j+1] - v[i,j-1])
		- c * dx * dy ;
