var x{1..3} integer >= 0 <= 4;

maximize o: x[1] + 3*x[2] + 6*x[3];

s.t. c1: alldiff(x[3], x[2], x[1]);
s.t. c2: abs(x[3] - x[3]) > 2;
s.t. c3: ((if x[1] >= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
