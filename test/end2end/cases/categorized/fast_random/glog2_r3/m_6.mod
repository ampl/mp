var x{1..3} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: (x[2] <= 3) ==> (x[1] >= 3) else (x[2] != 3);
s.t. c2: (x[3] < 1) or (x[3] < 0) or (x[3] < 1);
