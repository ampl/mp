var x{1..3} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: not ((x[3] <= 2) and (x[1] >= 2));
s.t. c2: abs(x[2] - x[1]) != 1;
s.t. c3: ((if x[1] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
