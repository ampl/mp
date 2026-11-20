var x{1..3} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (4)*x[2] + (3)*x[3];
s.t. c1: numberof 4 in (x[1], x[2], x[3]) >= 0;
s.t. c2: max(x[2], x[2]) - min(x[1], x[3]) < 3;
