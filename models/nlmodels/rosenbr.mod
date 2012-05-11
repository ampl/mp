# rosenbr.mod	OUR2-AN-2-0
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Rosenbrock function

# Ref.: J. J. More', B. S. Garbow, K. E. Hillstrom, "Testing Unconstrained
# Optimization Software", ACM Transactions on Mathematical Software,
# vol.7, no.1, 1981, pp. 17-41.

# Original reference:
# H. H. Rosenbrock, "An Automatic Method for Finding the Greatest or Least
# Value of a Function", Computer J., v. 3, 1960, pp. 175-184.

# Number of variables:  2
# Number of constraints:  0
# Objective nonseparable, nonconvex


var x{1..2};
var f1 = 10*(x[2] - x[1]^2);
var f2 = 1 - x[1];

minimize norm:
     f1^2 + f2^2;

data;
var x:=
        1  -1.2
        2   1;
