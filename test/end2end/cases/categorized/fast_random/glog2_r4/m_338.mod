var x{1..5} integer >= 0 <= 1;

minimize o: 2*x[1] + 6*x[2] + 3*x[3] + x[4] + 3*x[5];

s.t. c1: (x[2] == 1) ==> (x[4] >= 1) else (x[1] <= 1);
