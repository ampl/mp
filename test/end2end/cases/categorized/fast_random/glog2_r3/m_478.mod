var x{1..6} integer >= 0 <= 1;

maximize o: -x[1] - 5*x[2] - x[3] - 2*x[4] + 2*x[5] + 5*x[6];

s.t. c1: (if x[1] < 0 then x[5] else x[4]) != 0;
s.t. c2: (x[3] >= 1) or (x[5] > 0) or (x[2] >= 0);
s.t. c3: (if x[4] >= 0 then x[2] else x[3]) == 1;
s.t. c4: ((if x[1] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
