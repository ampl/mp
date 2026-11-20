var x{1..3} integer >= 0 <= 2;

minimize o: 0*x[1] - 3*x[2] + 2*x[3];

s.t. c1: (x[2] >= 1) ==> (x[1] <= 2) else (x[3] < 0);
s.t. c2: (x[2] > 1) or (x[2] != 1) or (x[3] <= 0);
