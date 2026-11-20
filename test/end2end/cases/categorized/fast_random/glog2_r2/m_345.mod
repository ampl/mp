var x{1..3} integer >= 0 <= 2;
minimize o: (6)*x[1] + (2)*x[2] + (-3)*x[3];
s.t. c1: numberof 1 in (x[1], x[2], x[3]) != 3;
s.t. c2: abs(x[3] - x[2]) > 1;
s.t. c3: max(x[3], x[2]) - min(x[1], x[1]) < 0;
s.t. c4: (x[1] > 1) <==> (x[3] + x[2] >= 2);
s.t. c5: (x[3] <= 0) or (x[2] < 1) or (x[2] != 0);
