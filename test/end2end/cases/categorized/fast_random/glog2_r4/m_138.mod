var x{1..6} integer >= 0 <= 1;

minimize o: -3*x[1] + x[2] + 3*x[3] + 5*x[4] - 5*x[5] + 0*x[6];

s.t. c1: (x[6] != 1) and (x[2] < 1) ==> (x[2] >= 1);
