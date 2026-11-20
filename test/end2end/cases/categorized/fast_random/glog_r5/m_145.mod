var x{1..5} integer >= -2 <= 2;

maximize o: -3*x[1] + x[2] - 3*x[3] + x[4] + 3*x[5];

s.t. c1: (x[2] < 1) or (x[3] != 2) or (x[5] > 1);
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) <= 1;
s.t. c3: abs(x[3] - x[1]) < 2;
s.t. c4: ((if x[4] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
