var x{1..3} integer >= -2 <= 1;
minimize o: (3)*x[1] + (-1)*x[2] + (-2)*x[3];
s.t. c1: numberof 0 in (x[1], x[2], x[3]) <= 0;
