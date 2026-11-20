var x{1..3} integer >= -2 <= 3;

minimize o: 0*x[1] - 4*x[2] + x[3];

s.t. c1: (x[1] < 0) ==> (x[3] <= 3);
s.t. c2: not ((x[3] != -1) and (x[2] < 2));
