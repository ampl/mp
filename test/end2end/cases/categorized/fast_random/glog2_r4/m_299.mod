var x{1..3} integer >= 0 <= 4;

maximize o: -3*x[1] - 3*x[2] + x[3];

s.t. c1: ((if x[2] >= 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c2: max(x[3], x[3]) - min(x[1], x[3]) != 4;
s.t. c3: ((if x[2] >= 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c4: max(x[3], x[3]) - min(x[1], x[3]) != 4;
