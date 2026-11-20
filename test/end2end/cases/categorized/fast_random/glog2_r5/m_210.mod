var x{1..4} integer >= 0 <= 4;

minimize o: 0*x[1] + 3*x[2] - 3*x[3] - 3*x[4];

s.t. c1: (x[2] > 1) ==> (x[4] == 1);
