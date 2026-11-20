var x{1..4} integer >= 0 <= 1;

minimize o: -4*x[1] - 2*x[2] + 5*x[3] + 5*x[4];

s.t. c1: (x[1] >= 1) and (x[1] >= 0) ==> (x[2] >= 1);
