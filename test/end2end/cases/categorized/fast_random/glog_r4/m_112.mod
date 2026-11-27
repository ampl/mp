var x{1..5} integer >= 0 <= 4;

minimize o: -x[1] - 2*x[2] - 4*x[3] + 5*x[4] - 3*x[5];

s.t. c1: abs(x[2] - x[2]) < 1;
s.t. c2: (x[2] < 0) and (x[4] < 0) ==> (x[1] >= 0);
