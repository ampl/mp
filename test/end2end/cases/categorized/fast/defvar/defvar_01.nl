g3 0 1 0	# problem defvar_01
 2 1 1 0 0	# vars, constraints, objectives, ranges, eqns
 0 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 2 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 2 2	# nonzeros in Jacobian, gradients
 2 1	# max name lengths: constraints, variables
 0 0 0 0 3	# common exprs: b,c,o,c1,o1
b	#2 bounds (on variables)
1 0
2 0
r	#1 ranges (rhs's)
2 3
C0	#L1
n0
V2 1 2	#dx
0 0.1
n2.4
V3 1 2	#dy
1 0.2
n8
V4 0 2	#dxy
o0	# + 
o2	#*
o1	# - 
n2.6
v2	#dx
o2	#*
n1.2
v3	#dy
n0.8
O0 0	#O1
v4	#dxy
k1	#intermediate Jacobian column lengths
1
J0 2
0 -1
1 1
G0 2
0 0
1 0
