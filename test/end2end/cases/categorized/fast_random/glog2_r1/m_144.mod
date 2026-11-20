var x{1..5} integer >= -2 <= 2;

minimize o: -5*x[1] - 4*x[2] + x[3] + 4*x[4] + x[5];

s.t. c1: (x[2] > -1) ==> (x[4] > 2) else (x[1] <= -1);
