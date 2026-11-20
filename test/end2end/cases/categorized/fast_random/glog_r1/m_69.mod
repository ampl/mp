var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] + 4*x[2] + 2*x[3] - 5*x[4];

s.t. c1: (x[1] >= 2) <==> (x[4] + x[3] > 2);
s.t. c2: alldiff(x[2], x[4], x[3], x[1]);
s.t. c3: ((if x[2] != 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
