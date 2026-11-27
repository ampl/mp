var x{1..5} integer >= 0 <= 2;

minimize o: -5*x[1] - x[2] + x[3] + x[4] + x[5];

s.t. c1: (if x[1] < 2 then x[5] else x[3]) > 0;
s.t. c2: max(x[3], x[5]) - min(x[3], x[5]) <= 2;
s.t. c3: not ((x[3] < 1) and (x[5] <= 2));
s.t. c4: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
