g3 2 1 0	# problem viold_with_obj
 5 2 1 0 0 1	# vars, constraints, objectives, ranges, eqns, lcons
 2 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 5 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 2 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 5 0	# nonzeros in Jacobian, gradients
 47 28	# max name lengths: constraints, variables
 0 0 0 3 0	# common exprs: b,c,o,c1,o1
b	#5 bounds (on variables)
0 0 229
0 109.20000000000002 364.00000000000006
0 109.20000000000002 364.00000000000006
0 0 1
0 0 1
r	#2 ranges (rhs's)
1 1500
1 0
C0	#charge_diff_const['2024-10-16 21:00:00']
v6	#charge_diff_var['2024-10-16 21:00:00']
C1	#charge_diff_const['2024-10-16 21:15:00']
o1	# - 
v7	#charge_diff_var['2024-10-16 21:15:00']
o35	# if 
o21	# && 
o24	# == 
v4	#y_ind['2024-10-16 21:15:00']
n1
o24	# == 
v3	#y_ind['2024-10-16 21:00:00']
n1
n50
n1500
V5 1 3	#cooling_residual['2024-10-16 20:45:00']
0 -1
n229
L0	#CustLimitationConstraint['2024-10-16 20:45:00']
o20	# || 
o20	# || 
o34	#!
o29	# > 
v5	#cooling_residual['2024-10-16 20:45:00']
n0
o34	#!
o28	# >= 
v0	#x['2024-10-16 20:45:00']
n15
o28	# >= 
v5	#cooling_residual['2024-10-16 20:45:00']
n70
O0 1	#Obj
n0
k4	#intermediate Jacobian column lengths
0
2
3
4
J0 1
1 0
J1 4
1 0
2 0
3 0
4 0
G0 4
1 1
2 -1
3 1000
4 1000
