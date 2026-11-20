var x{1..5} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-4)*x[2] + (3)*x[3] + (0)*x[4] + (1)*x[5];
s.t. c1: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c2: (if x[3] <= 3 then x[2] else x[4]) >= 3;
s.t. c3: (x[4] <= 1) or (x[4] > 2) or (x[5] < 2);
s.t. c4: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c5: (x[4] <= 1) or (x[4] > 2) or (x[5] < 2);
