var x{1..4} integer >= -2 <= 1;

minimize o: -x[1] - 2*x[2] - 4*x[3] - x[4];

s.t. c1: (x[3] == 0) ==> (x[2] >= 0);
