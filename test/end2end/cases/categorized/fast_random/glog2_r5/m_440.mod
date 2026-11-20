var x{1..5} integer >= 0 <= 4;

minimize o: 0*x[1] - 3*x[2] - 2*x[3] - 5*x[4] + 2*x[5];

s.t. c1: not ((x[2] < 1) and (x[5] <= 2));
s.t. c2: ((if x[2] < 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
