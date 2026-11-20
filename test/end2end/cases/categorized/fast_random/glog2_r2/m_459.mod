var x{1..5} integer >= 0 <= 2;

minimize o: 5*x[1] + 2*x[2] + 4*x[3] + 6*x[4] - x[5];

s.t. c1: alldiff(x[2], x[5], x[1], x[4]);
s.t. c2: (x[3] != 1) <==> (x[3] + x[5] >= 1);
s.t. c3: max(x[2], x[1]) - min(x[5], x[1]) > 0;
s.t. c4: (x[2] != 2) and (x[5] <= 1) ==> (x[1] >= 2);
s.t. c5: alldiff(x[2], x[5], x[1], x[4]);
s.t. c6: (x[2] != 2) and (x[5] <= 1) ==> (x[1] >= 2);
