var x{1..5} integer >= 0 <= 1;

maximize o: 2*x[1] + 4*x[2] - 5*x[3] - 3*x[4] - 2*x[5];

s.t. c1: max(x[2], x[5]) - min(x[2], x[1]) <= 0;
s.t. c2: (x[1] <= 0) <==> (x[5] + x[3] <= 1);
s.t. c3: ((if x[3] < 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
