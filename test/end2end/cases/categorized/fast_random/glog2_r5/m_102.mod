var x{1..3} integer >= -2 <= 0;
minimize o: (-5)*x[1] + (-4)*x[2] + (2)*x[3];
s.t. c1: numberof -2 in (x[1], x[2], x[3]) < 1;
