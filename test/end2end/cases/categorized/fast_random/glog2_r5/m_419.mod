var x{1..6} integer >= -2 <= 2;

maximize o: -2*x[1] + 6*x[2] - 5*x[3] + 2*x[4] + 2*x[5] + 6*x[6];

s.t. c1: (if x[2] >= 2 then x[3] else x[4]) == 1;
s.t. c2: ((if x[4] > 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: (if x[2] >= 2 then x[3] else x[4]) == 1;
s.t. c4: ((if x[4] > 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
