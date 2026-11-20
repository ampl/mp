var x{1..6} integer >= 0 <= 2;

minimize o: -5*x[1] - 2*x[2] + 5*x[3] - 4*x[4] + 3*x[5] + 3*x[6];

s.t. c1: max(x[2], x[6]) - min(x[4], x[3]) > 1;
s.t. c2: (x[1] <= 0) and (x[2] == 2) ==> (x[2] >= 0);
