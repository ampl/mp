g3 1 1 0	# problem ch3
 3 3 1 0 3	# vars, constraints, objectives, ranges, eqns
 2 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 3 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 9 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 9 0 0 0 3	# common exprs: b,c,o,c1,o1
V3 1 0	#T[1,1]
0 2
n-1
V4 1 0	#T[1,2]
1 2
n-1
V5 1 0	#T[1,3]
2 2
n-1
V6 0 0	#T[2,1]
o0	# + 
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v0	#x[1]
n1
v3	#T[1,1]
n-1
V7 0 0	#T[2,2]
o0	# + 
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v1	#x[2]
n1
v4	#T[1,2]
n-1
V8 0 0	#T[2,3]
o0	# + 
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v2	#x[3]
n1
v5	#T[1,3]
n-1
V9 0 0	#nl(T[3,1])
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v0	#x[1]
n1
v6	#T[2,1]
V10 0 0	#nl(T[3,2])
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v1	#x[2]
n1
v7	#T[2,2]
V11 0 0	#nl(T[3,3])
o2	#*
o2	#*
n2
o1	# - 
o2	#*
n2
v2	#x[3]
n1
v8	#T[2,3]
C0	#eqn[2]
o54	#sumlist
3
o2	#*
n0.3333333333333333
v8	#T[2,3]
o2	#*
n0.3333333333333333
v7	#T[2,2]
o2	#*
n0.3333333333333333
v6	#T[2,1]
C1	#eqn[3]
o54	#sumlist
3
o2	#*
n0.3333333333333333
v11	#nl(T[3,3])
o2	#*
n0.3333333333333333
v10	#nl(T[3,2])
o2	#*
n0.3333333333333333
v9	#nl(T[3,1])
C2	#eqn[1]
n0
V12 1 4	#T[3,1]
0 -2
o0	# + 
v9	#nl(T[3,1])
n1
V13 1 4	#T[3,2]
1 -2
o0	# + 
v10	#nl(T[3,2])
n1
V14 1 4	#T[3,3]
2 -2
o0	# + 
v11	#nl(T[3,3])
n1
O0 0	#ssq
o54	#sumlist
3
o2	#*
n0.5
o5	#^
o2	#*
n0.3333333333333333
o54	#sumlist
3
v3	#T[1,1]
v4	#T[1,2]
v5	#T[1,3]
n2
o2	#*
n0.5
o5	#^
o1	# - 
o2	#*
n0.3333333333333333
o54	#sumlist
3
v6	#T[2,1]
v7	#T[2,2]
v8	#T[2,3]
n-0.3333333333333333
n2
o2	#*
n0.5
o5	#^
o2	#*
n0.3333333333333333
o54	#sumlist
3
v12	#T[3,1]
v13	#T[3,2]
v14	#T[3,3]
n2
x3	# initial guess
0 0.25
1 0.5
2 0.75
r	#3 ranges (rhs's)
4 -0.3333333333333333
4 -1
4 1
b	#3 bounds (on variables)
3
3
3
k2	#intermediate Jacobian column lengths
3
6
J0 3
0 0
1 0
2 0
J1 3
0 -0.6666666666666666
1 -0.6666666666666666
2 -0.6666666666666666
J2 3
0 0.6666666666666666
1 0.6666666666666666
2 0.6666666666666666
G0 3
0 0
1 0
2 0
