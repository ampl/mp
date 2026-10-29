g3 1 1 0	# problem inline_01
 7 3 1 0 0	# vars, constraints, objectives, ranges, eqns
 3 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 6 7 5	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 5 1 1	# discrete variables: binary, integer, nonlinear (b,c,o)
 10 7	# nonzeros in Jacobian, gradients
 12 4	# max name lengths: constraints, variables
 1 1 0 0 0	# common exprs: b,c,o,c1,o1
b	#7 bounds (on variables)
0 0 1
0 0 1
0 0 4
0 0 4
0 0 4
0 0 4
0 0 1
r	#3 ranges (rhs's)
1 5
1 0
2 1
V7 0 0	#count2
o59	#count
3
o24	# == 
v2	#x[1]
n1
o24	# == 
v3	#x[2]
n2
o24	# == 
v4	#x[3]
n3
V8 0 0	#if2
o35	# if 
o30	# != 
v1	#b[3]
n0
n2
n0
C0	#C1
o54	#sumlist
3
v8	#if2
o16	#-
o60	#numberof
4
v5	#x[4]
v2	#x[1]
v3	#x[2]
v4	#x[3]
o59	#count
2
o34	#!
o30	# != 
v0	#b[2]
n0
o34	#!
o30	# != 
v1	#b[3]
n0
C1	#Fix_If2
v8	#if2
C2	#Bound_Count2
v7	#count2
O0 0	#Obj
o54	#sumlist
3
o16	#-
v7	#count2
o35	# if 
o30	# != 
v6	#b[1]
n0
n5
n2
o59	#count
2
o34	#!
o30	# != 
v0	#b[2]
n0
o34	#!
o30	# != 
v1	#b[3]
n0
k6	#intermediate Jacobian column lengths
1
3
5
7
9
10
J0 6
0 0
1 0
2 0
3 0
4 0
5 0
J1 1
1 0
J2 3
2 0
3 0
4 0
G0 7
0 0
1 0
2 1
3 -1
4 1
5 -1
6 0
