var x{1..3} integer >= -2 <= 3;
minimize o: (-4)*x[1] + (-5)*x[2] + (0)*x[3];
s.t. c1: numberof 2 in (x[1], x[2], x[3]) < 1;
s.t. c2: not ((x[3] != 1) and (x[3] == 0));
