var x{1..6} integer >= 0 <= 2;

minimize o: 0*x[1] + 6*x[2] - 4*x[3] - 4*x[4] + 5*x[5] - 2*x[6];

s.t. c1: (x[6] < 0) or (x[5] > 0) or (x[5] != 0);
s.t. c2: (x[5] < 0) ==> (x[6] > 1) else (x[2] < 2);
