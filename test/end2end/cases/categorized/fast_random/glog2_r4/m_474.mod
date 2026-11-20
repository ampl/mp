var x{1..5} integer >= -3 <= 2;

maximize o: 4*x[1] - x[2] + 0*x[3] + 6*x[4] + 6*x[5];

s.t. c1: (x[1] > 1) ==> (x[5] >= -1) else (x[2] > 0);
s.t. c2: alldiff(x[2], x[3], x[1], x[4]);
s.t. c3: not ((x[4] <= 0) and (x[5] > 0));
s.t. c4: (x[5] <= -2) <==> (x[3] + x[4] > 1);
