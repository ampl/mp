var x{1..4} integer >= 0 <= 4;

minimize o: -4*x[1] - 5*x[2] + 3*x[3] + 5*x[4];

s.t. c1: (x[3] != 2) and (x[1] <= 2) ==> (x[1] >= 2);
