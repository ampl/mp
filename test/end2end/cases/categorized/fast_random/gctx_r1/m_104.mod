var x{1..5} integer >= 0 <= 2;

minimize o: -x[1] - x[2] + x[3] + 0*x[4] - 2*x[5];

s.t. c1: (x[2] <= 0) ==> (x[2] <= 0) else (x[1] < 0);
s.t. c2: (x[2] <= 0) or (x[2] <= 0);
