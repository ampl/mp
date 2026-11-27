var x{1..4} integer >= -2 <= 1;

maximize o: -5*x[1] - 3*x[2] - 2*x[3] - 3*x[4];

s.t. c1: ((if x[4] <= -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[2], x[3], x[4], x[1]);
s.t. c3: ((if x[4] <= -2 then 1 else 0) + (if x[4] >= -2 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[2], x[3], x[4], x[1]);
