var x{1..5} integer >= 0 <= 4;

minimize o: -4*x[1] + 4*x[2] + x[3] - x[4] + x[5];

s.t. c1: (x[3] > 2) ==> (x[1] <= 0);
