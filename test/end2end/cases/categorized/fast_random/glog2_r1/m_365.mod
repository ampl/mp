var x{1..6} integer >= 0 <= 2;

minimize o: 2*x[1] - 2*x[2] - 3*x[3] + 0*x[4] + x[5] - 3*x[6];

s.t. c1: (x[5] >= 1) ==> (x[3] <= 1);
s.t. c2: (x[1] > 2) ==> (x[5] <= 2) else (x[1] != 2);
s.t. c3: alldiff(x[5], x[3], x[4]);
s.t. c4: abs(x[1] - x[1]) > 2;
s.t. c5: alldiff(x[4], x[2], x[1]);
s.t. c6: ((if x[6] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
