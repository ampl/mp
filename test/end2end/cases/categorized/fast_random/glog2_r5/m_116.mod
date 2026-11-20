var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: max(x[1], x[1]) - min(x[2], x[3]) <= 0;
s.t. c2: ((if x[2] < 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c3: (x[4] > 2) <==> (x[2] + x[4] < 3);
s.t. c4: alldiff(x[3], x[4]);
