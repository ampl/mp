var x{1..6} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (-5)*x[2] + (0)*x[3] + (-2)*x[4] + (-4)*x[5] + (2)*x[6];
s.t. c1: not ((x[4] == 0) and (x[6] < 2));
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 0;
s.t. c3: (if x[3] < 1 then x[3] else x[2]) >= 2;
