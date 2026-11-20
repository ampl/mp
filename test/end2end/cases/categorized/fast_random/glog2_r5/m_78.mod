var x{1..5} integer >= 0 <= 1;

minimize o: 6*x[1] + 3*x[2] + 4*x[3] + 5*x[4] + 4*x[5];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
s.t. c2: (if x[5] == 0 then x[3] else x[3]) <= 1;
s.t. c3: abs(x[1] - x[5]) <= 1;
s.t. c4: ((if x[5] < 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
