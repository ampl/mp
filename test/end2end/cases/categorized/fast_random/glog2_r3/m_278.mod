var x{1..5} integer >= 0 <= 1;

maximize o: -5*x[1] + 0*x[2] - x[3] - 2*x[4] + 0*x[5];

s.t. c1: alldiff(x[3], x[1], x[5], x[4]);
s.t. c2: ((if x[3] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..5} x[i] == 1;
