var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] - 5*x[2] + 0*x[3] - 3*x[4];

s.t. c1: (x[1] != 1) ==> (x[1] <= 3);
