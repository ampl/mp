var x{1..4} integer >= -2 <= 0;

minimize o: -x[1] - 5*x[2] + 0*x[3] + x[4];

s.t. c1: (x[3] == -1) ==> (x[1] >= 0);
