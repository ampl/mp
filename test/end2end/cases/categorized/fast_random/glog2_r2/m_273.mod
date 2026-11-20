var x{1..6} integer >= 0 <= 4;
minimize o: (-1)*x[1] + (0)*x[2] + (-1)*x[3] + (-4)*x[4] + (6)*x[5] + (1)*x[6];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
s.t. c2: max(x[5], x[1]) - min(x[1], x[2]) < 1;
s.t. c3: ((if x[5] == 4 then 1 else 0) + (if x[5] >= 4 then 1 else 0)) >= 1;
s.t. c4: ((if x[6] >= 4 then 1 else 0) + (if x[6] >= 4 then 1 else 0)) >= 1;
s.t. c5: x[2] == 4;
s.t. c6: x[3] == 4;
