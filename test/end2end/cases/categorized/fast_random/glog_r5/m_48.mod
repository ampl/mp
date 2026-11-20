var x{1..4} integer >= 0 <= 2;

minimize o: -2*x[1] + 0*x[2] - x[3] + x[4];

s.t. c1: (x[4] < 2) ==> (x[4] > 0);
