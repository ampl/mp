var x{1..4} integer >= 0 <= 1;

minimize o: 2*x[1] + 5*x[2] + 3*x[3] - x[4];

s.t. c1: (x[1] < 0) ==> (x[2] > 0) else (x[3] <= 1);
