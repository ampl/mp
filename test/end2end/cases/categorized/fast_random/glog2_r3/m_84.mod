var x{1..3} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (-3)*x[2] + (6)*x[3];
s.t. c1: numberof -2 in (x[1], x[2], x[3]) == 1;
s.t. c2: (x[1] != -1) ==> (x[1] < 2);
