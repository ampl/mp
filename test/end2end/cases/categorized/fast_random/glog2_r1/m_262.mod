var x{1..4} integer >= 0 <= 2;

minimize o: -5*x[1] + 6*x[2] - 4*x[3] - 4*x[4];

s.t. c1: (if x[1] >= 1 then x[4] else x[3]) >= 0;
s.t. c2: ((if x[2] < 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..4} x[i] == 2;
s.t. c4: (if x[2] == 1 then x[4] else x[4]) == 1;
