var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] - x[2] - x[3] - 2*x[4];

s.t. c1: (x[4] > 0) ==> (x[2] <= 2) else (x[4] == 0);
