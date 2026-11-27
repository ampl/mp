var x{1..7} integer >= 0 <= 1;

minimize o: -5*x[1] - 2*x[2] + 6*x[3] + 6*x[4] + 0*x[5] - 2*x[6] - 4*x[7];

s.t. c1: (x[6] != 0) ==> (x[4] >= 1);
s.t. c2: (x[4] > 0) and (x[5] <= 1) ==> (x[4] >= 0);
s.t. c3: alldiff(x[6], x[1], x[4], x[7]);
s.t. c4: alldiff(x[2], x[4], x[7], x[5]);
