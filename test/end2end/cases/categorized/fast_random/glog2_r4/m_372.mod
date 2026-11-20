var x{1..4} integer >= -2 <= 2;

minimize o: 4*x[1] + 4*x[2] + 4*x[3] + x[4];

s.t. c1: (x[1] > -2) ==> (x[1] <= 1);
s.t. c2: not ((x[3] >= -2) and (x[3] > 0));
