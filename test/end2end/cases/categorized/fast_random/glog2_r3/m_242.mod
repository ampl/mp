var x{1..6} integer >= -2 <= 1;

minimize o: 3*x[1] + x[2] + 2*x[3] - 2*x[4] + 6*x[5] + 6*x[6];

s.t. c1: (x[2] > -1) or (x[1] > 1) or (x[1] < 1);
s.t. c2: max(x[6], x[3]) - min(x[6], x[6]) <= -2;
s.t. c3: ((if x[5] > 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c4: (if x[2] >= 1 then x[3] else x[4]) != 1;
