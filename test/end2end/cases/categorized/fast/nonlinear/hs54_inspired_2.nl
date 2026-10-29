g3 2 1 0	# problem hs54_inspired_2
 6 1 1 0 1	# vars, constraints, objectives, ranges, eqns
 0 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 0 6 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 2 6	# nonzeros in Jacobian, gradients
 3 4	# max name lengths: constraints, variables
 0 0 0 0 1	# common exprs: b,c,o,c1,o1
b	#6 bounds (on variables)
0 0 200
0 -10 10
0 -5 12
0 -15 20
0 -100 1
0 -5 22
r	#1 ranges (rhs's)
4 176
C0	#C1
n0
V6 0 2	#h
o54	#sumlist
5
o2	#*
n0.00021258503401360543
o2	#*
o54	#sumlist
3
o3	#/
o5	#^
o0	# + 
n-10
v0	#x[1]
n2
n6400
o3	#/
o2	#*
o0	# + 
n-10
v0	#x[1]
o0	# + 
n-1
v1	#x[2]
n200
o3	#/
o5	#^
o0	# + 
n-1
v1	#x[2]
n2
n3000
o5	#^
o0	# + 
n-200
v2	#x[3]
n2
o2	#*
n0.004
o5	#^
o0	# + 
n-10
v3	#x[4]
n2
o2	#*
n0.004
o5	#^
o0	# + 
n-0.1
v4	#x[5]
n2
o2	#*
n0.00021258503401360543
o2	#*
o2	#*
o54	#sumlist
4
o3	#/
o2	#*
o0	# + 
n-20
v0	#x[1]
v0	#x[1]
n6400
o3	#/
o2	#*
o0	# + 
n-10
v0	#x[1]
v1	#x[2]
n200
o3	#/
o2	#*
o54	#sumlist
3
o3	#/
o5	#^
o0	# + 
n-10
v2	#x[3]
n2
n6400
o3	#/
o2	#*
o0	# + 
n-10
v1	#x[2]
o0	# + 
n-1
v1	#x[2]
n200
o3	#/
o5	#^
o0	# + 
n-1
v3	#x[4]
n2
n3000
o5	#^
o0	# + 
n-200
v2	#x[3]
n2
n4704
o3	#/
o2	#*
o0	# + 
n-1
v1	#x[2]
v1	#x[2]
n3000
o0	# + 
n-200
v2	#x[3]
v3	#x[4]
o2	#*
n0.0004
o5	#^
o0	# + 
n-100
v5	#x[6]
n2
O0 0	#Obj
o16	#-
o44	#exp
o16	#-
o3	#/
v6	#h
n2
k5	#intermediate Jacobian column lengths
1
2
2
2
2
J0 2
0 1
1 40
G0 6
0 0
1 0
2 0
3 0
4 0
5 0
