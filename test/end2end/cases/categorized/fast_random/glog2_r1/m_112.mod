var x{1..6} integer >= 0 <= 3;

minimize o: -x[1] - 5*x[2] - x[3] + 6*x[4] - 4*x[5] - 5*x[6];

s.t. c1: (x[4] <= 2) or (x[6] > 1) or (x[1] == 3);
s.t. c2: not ((x[4] != 0) and (x[4] <= 0));
s.t. c3: ((if x[5] < 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
