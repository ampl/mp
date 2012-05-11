# bertsekas.mod

# This problem is based on a traffic assigment problem, found in:
# D.P. Bertsekas and E.M. Gafni, "Projection Methods for Variational
# Inequalities with Application to the Traffic Assignment Problem",
# Mathematical Programming Study 17(1982), 139-159.

# Just so you know, x(i) is the flow (from node(i) to node(i+3)) on the
# (longer) outside loop, z(i) the flow on the inside loop, and y the
# dual variables.

# The source problem here is a VI, where the constraints consist of
# non-negative bounds and x(I) + z(I) = demand(I).  Without much
# fiddling, we can write the problem as a MCP.

# We can write the problem in (at least) 2 different ways.  In the
# first, y is free, and the demand constraints are satisfied exactly.
# In the second, we change the demand constraint to an inequality, where
# extra flow is allowed.  This doesn't occur in solutions, and we can
# write the problem as an NCP.

set Nodes  circular;	# nodes in the network
set Types;		# types of arcs

param demand {Nodes};	# flow demand, i -> i+3
param c {Nodes,Types};	# delay coefficients for links in the network

param gamma;	# gamma is a coupling coefficient.  The delay on some of
		# the links is increased  gamma * delay on merging links.

var x {Nodes};	# outside flow
var z {Nodes};	# inside flow
var y {Nodes};	# dual vars

var fx {i in Nodes} = 1 + x[i] + x[i]^2;
var fz {i in Nodes} = 1 + z[i] + z[i]^2;

s.t. delay_o {i in Nodes}: 0 <= x[i]
     complements
	  c[i,"oon"]*fx[i]
	+ gamma * (1 + x[next(i,Nodes,3)] + x[next(i,Nodes,4)]
		   + (x[next(i,Nodes,3)]+x[next(i,Nodes,4)])^2)
	+ 10 * c[i,"ohy"]*(1 + x[i] + x[next(i,Nodes,3)]
			   + x[next(i,Nodes,4)]
			   + (x[i] + x[next(i,Nodes,3)]
				+ x[next(i,Nodes,4)])^2)
	+ 2 * gamma * fx[next(i,Nodes,3)]
	+ c[next(i),"oby"] * (1 + x[i] + x[next(i,Nodes,4)]
			      + (x[i] + x[next(i,Nodes,4)])^2)
	+ 10 * c[next(i),"ohy"]
	     * (1 + x[i] + x[next(i)] + x[next(i,Nodes,4)]
		+ (x[i] + x[next(i)] + x[next(i,Nodes,4)])^2)
	+ 2 * gamma * fx[next(i,Nodes,4)]
	+ c[next(i,Nodes,2),"oby"]
	  * (1 + x[i] + x[next(i)] + (x[i]+x[next(i)])^2)
	+ 10 * c[next(i,Nodes,2),"ohy"]
	     * (1 + x[i] + x[next(i)] + x[next(i,Nodes,2)]
		  + (x[i] + x[next(i)] + x[next(i,Nodes,2)])^2)
	+ 2 * gamma * fx[i]
	+ c[next(i,Nodes,3),"ooff"] * fx[i]
	>= y[i];


s.t. delay_i {i in Nodes}: 0 <= z[i]
     complements
	  c[i,"ion"]*fz[i]
	+ gamma * fz[next(i)]
	+ 10 * c[i,"ihy"]*(1 + z[i] + z[next(i)] + (z[i] + z[next(i)])^2)
	+ 2 * gamma * fz[next(i)]
	+ c[prev(i),"iby"]*fz[i]
	+ 10 * c[prev(i),"ihy"] * (1 + z[i] + z[next(i,Nodes,4)]
				   + (z[i]+z[next(i,Nodes,4)])^2)
	+ 2 * gamma * fz[i]
	+ c[prev(i,Nodes,2),"ioff"] * fz[i]
	>= y[i];

s.t. dem_cons {i in Nodes}: x[i] + z[i] == demand[i];
