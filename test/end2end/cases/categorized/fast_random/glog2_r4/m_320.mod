var x{1..6} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-2)*x[2] + (-5)*x[3] + (0)*x[4] + (-1)*x[5] + (-5)*x[6];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
s.t. c2: (if x[2] < 1 then x[6] else x[4]) > 1;
s.t. c3: alldiff(x[6], x[2], x[5]);
