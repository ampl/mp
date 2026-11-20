var x{1..5} integer >= 0 <= 2;

minimize o: 0*x[1] + 2*x[2] + 3*x[3] + 6*x[4] - 3*x[5];

s.t. c1: alldiff(x[5], x[3]);
s.t. c2: (x[4] == 1) and (x[4] == 1) ==> (x[5] >= 1);
s.t. c3: (x[2] > 1) or (x[5] <= 0) or (x[5] >= 1);
s.t. c4: alldiff(x[5], x[3]);
s.t. c5: (x[2] > 1) or (x[5] <= 0) or (x[5] >= 1);
