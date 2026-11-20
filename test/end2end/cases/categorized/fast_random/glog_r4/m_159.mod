var x{1..3} integer >= 0 <= 2;

minimize o: 5*x[1] - 3*x[2] - 4*x[3];

s.t. c1: (x[2] != 0) and (x[1] < 1) ==> (x[3] >= 0);
