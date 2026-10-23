g3 2 1 0	# problem power5
 2 3 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 2 2	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 6 2	# nonzeros in Jacobian, gradients
 3 4	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b	#2 bounds (on variables)
2 2.8
2 2.8
r	#3 ranges (rhs's)
2 1
2 0.5
1 10.2
C0	#c2
o0	# + 
o2	#*
n0.1111111111111111
o5	#^
o0	# + 
n-3.4
o0	# + 
o2	#*
n2.3
v0	#x[1]
o2	#*
n1.57
v1	#x[2]
n5
o5	#^
v1	#x[2]
n2
C1	#c3
o2	#*
n8
o5	#^
v0	#x[1]
n2
C2	#c1
n0
O0 1	#Sum
o0	# + 
o2	#*
n-5
o5	#^
o0	# + 
n-0.7
v0	#x[1]
n2
o5	#^
v1	#x[2]
n7
k1	#intermediate Jacobian column lengths
3
J0 2
0 0
1 0
J1 2
0 0
1 1
J2 2
0 2
1 1
G0 2
0 0
1 0
