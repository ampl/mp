var x{1..6} integer >= 0 <= 2;

minimize o: 5*x[1] - 3*x[2] + x[3] + 2*x[4] + 6*x[5] + 4*x[6];

s.t. c1: (x[6] <= 1) ==> (x[1] < 0);
