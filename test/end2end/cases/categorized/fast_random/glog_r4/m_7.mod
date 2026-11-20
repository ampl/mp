var x{1..3} integer >= 0 <= 2;

minimize o: -3*x[1] + 0*x[2] - 2*x[3];

s.t. c1: alldiff(x[3], x[2]);
s.t. c2: ((if x[3] == 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[2] - x[2]) <= 1;
s.t. c4: not ((x[3] == 2) and (x[2] >= 1));
