var x{1..4} integer >= 0 <= 1;

minimize o: -2*x[1] + x[2] - 4*x[3] + 5*x[4];

s.t. c1: (x[4] < 0) ==> (x[2] != 0);
