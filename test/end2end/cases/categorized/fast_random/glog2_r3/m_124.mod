var x{1..3} integer >= -3 <= 2;
minimize o: (-5)*x[1] + (4)*x[2] + (-5)*x[3];
s.t. c1: numberof -1 in (x[1], x[2], x[3]) >= 2;
