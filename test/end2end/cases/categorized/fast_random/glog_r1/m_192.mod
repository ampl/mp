var x{1..5} integer >= 0 <= 3;

minimize o: 2*x[1] + 2*x[2] - 3*x[3] + 0*x[4] + 0*x[5];

s.t. c1: abs(x[2] - x[4]) > 1;
s.t. c2: (x[1] <= 1) and (x[5] <= 0) ==> (x[4] >= 1);
