var x{1..6} integer >= 0 <= 4;

minimize o: 34*x[1] + 3000*x[2] - 1000*x[3] + 2000*x[4] + 51*x[5] - 999*x[6];

s.t. c1: not ((x[2] <= 4) and (x[1] == 2));
s.t. c2: ((if x[3] < 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
