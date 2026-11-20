var x{1..6} integer >= 0 <= 4;

minimize o: -x[1] + 6*x[2] - 5*x[3] + 2*x[4] - 5*x[5] + 6*x[6];

s.t. c1: (x[1] == 4) and (x[3] > 2) ==> (x[3] >= 4);
