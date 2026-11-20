var x{1..4} integer >= 0 <= 1;
minimize o: (-2)*x[1] + (-3)*x[2] + (5)*x[3] + (-3)*x[4];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c3: max(x[2], x[2]) - min(x[1], x[3]) == 1;
