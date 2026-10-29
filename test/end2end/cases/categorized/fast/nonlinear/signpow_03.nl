g3 0 1 0	# problem signpow_03
 3 1 1 0 0	# vars, constraints, objectives, ranges, eqns
 1 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 1 3 1	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 1 3	# nonzeros in Jacobian, gradients
 11 1	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b	#3 bounds (on variables)
0 -1.5 0.2
0 -2 3
0 -1.5 0.0025
r	#1 ranges (rhs's)
2 15
C0	#QuickInfeas
o0	# + 
o39	#sqrt
o2	#*
v0	#y
v0	#y
o39	#sqrt
o2	#*
o2	#*
v0	#y
v0	#y
v0	#y
O0 1	#Obj
o16	#-
o54	#sumlist
8
o2	#*
o2	#*
o15	#abs
v1	#x
v1	#x
v1	#x
o2	#*
o2	#*
o5	#^
o15	#abs
v0	#y
n2
o5	#^
o15	#abs
v0	#y
n0.3
o5	#^
v0	#y
n0.3
o16	#-
o2	#*
o2	#*
n360
v2	#z
v0	#y
o2	#*
o2	#*
o2	#*
o2	#*
o2	#*
o2	#*
o15	#abs
o5	#^
v0	#y
n0.22
o5	#^
v0	#y
n3
o5	#^
v0	#y
n-2
o15	#abs
v2	#z
o5	#^
o15	#abs
v2	#z
n4
v2	#z
o5	#^
v2	#z
n-4
o2	#*
o5	#^
o39	#sqrt
o5	#^
v1	#x
n2
n2
v1	#x
o2	#*
o2	#*
o5	#^
o39	#sqrt
o2	#*
v0	#y
v0	#y
n2
o5	#^
o39	#sqrt
o2	#*
o2	#*
v0	#y
v0	#y
v0	#y
n0.3
o5	#^
v0	#y
n3
o16	#-
o2	#*
o2	#*
n360
v2	#z
v0	#y
o2	#*
o2	#*
o2	#*
o2	#*
o2	#*
o2	#*
o15	#abs
o5	#^
v0	#y
n0.33
o5	#^
v0	#y
n3
o39	#sqrt
o5	#^
v0	#y
n4.4
o15	#abs
v2	#z
o5	#^
o39	#sqrt
o5	#^
v2	#z
n4
n2
v2	#z
o5	#^
v2	#z
n-4
k2	#intermediate Jacobian column lengths
1
1
J0 1
0 0
G0 3
0 26
1 4
2 0
