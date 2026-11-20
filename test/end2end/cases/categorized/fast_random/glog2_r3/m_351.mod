var x{1..3} integer >= -2 <= 0;
maximize o: (-5)*x[1] + (-2)*x[2] + (0)*x[3];
s.t. c1: abs(x[1] - x[2]) > 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3]) >= 1;
