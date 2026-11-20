var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] + 5*x[2] - 4*x[3] + 0*x[4];

s.t. c1: (x[3] == 1) and (x[1] > 1) ==> (x[2] >= 1);
s.t. c2: alldiff(x[2], x[3], x[4]);
