var x{1..6} integer >= -2 <= 2;

maximize o: -x[1] - 5*x[2] - 3*x[3] + 6*x[4] + 5*x[5] + 6*x[6];

s.t. c1: ((if x[2] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[5] <= -1) or (x[1] != 0) or (x[4] < 1);
s.t. c3: ((if x[1] != -1 then 1 else 0) + (if x[6] >= -1 then 1 else 0)) >= 1;
