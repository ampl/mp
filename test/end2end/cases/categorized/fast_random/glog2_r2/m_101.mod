var x{1..6} integer >= 0 <= 4;

minimize o: -4*x[1] - 2*x[2] + x[3] - x[4] + 6*x[5] - x[6];

s.t. c1: alldiff(x[4], x[5], x[1]);
s.t. c2: (x[6] < 5) and (x[1] > 0) ==> (x[5] >= 5);
