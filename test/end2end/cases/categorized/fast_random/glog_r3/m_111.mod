var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 2*x[2] + 5*x[3] + 6*x[4];

s.t. c1: max(x[1], x[4]) - min(x[3], x[4]) >= 1;
s.t. c2: (x[4] >= 1) ==> (x[4] > 4);
s.t. c3: ((if x[3] != 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
