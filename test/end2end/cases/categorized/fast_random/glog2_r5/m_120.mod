var x{1..6} integer >= 0 <= 4;

maximize o: -x[1] - x[2] - 5*x[3] - 3*x[4] + x[5] + 0*x[6];

s.t. c1: (x[4] <= 1) ==> (x[2] > 3) else (x[2] <= 2);
s.t. c2: alldiff(x[6], x[5], x[2], x[4]);
