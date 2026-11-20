var x{1..6} integer >= 0 <= 4;
minimize o: (-3)*x[1] + (-5)*x[2] + (0)*x[3] + (5)*x[4] + (-4)*x[5] + (2)*x[6];
s.t. c1: numberof 4 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
s.t. c2: not ((x[2] < 1) and (x[2] == 3));
s.t. c3: max(x[3], x[5]) - min(x[2], x[4]) >= 1;
