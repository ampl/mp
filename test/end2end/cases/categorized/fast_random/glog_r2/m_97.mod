var x{1..4} integer >= 0 <= 2;

minimize o: -2*x[1] + 0*x[2] - 3*x[3] + 3*x[4];

s.t. c1: not ((x[2] <= 2) and (x[4] >= 1));
s.t. c2: not ((x[4] != 2) and (x[2] != 0));
s.t. c3: (x[3] > 2) ==> (x[1] >= 2) else (x[4] < 0);
