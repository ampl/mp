var x{1..3} integer >= -3 <= 1;

minimize o: 4*x[1] + 6*x[2] - 2*x[3];

s.t. c1: (x[1] >= 0) and (x[1] > -1) ==> (x[2] >= 0);
