var x{1..4} integer >= 0 <= 1;

minimize o: 3*x[1] + 0*x[2] + 6*x[3] + x[4];

s.t. c1: ((if x[2] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[3], x[1]);
