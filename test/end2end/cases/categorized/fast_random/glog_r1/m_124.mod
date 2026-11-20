var x{1..4} integer >= -2 <= 0;
minimize o: (5)*x[1] + (-4)*x[2] + (0)*x[3] + (3)*x[4];
s.t. c1: (x[4] != 0) and (x[2] != -1) ==> (x[2] >= 0);
s.t. c2: (x[2] != -2) or (x[4] <= -2) or (x[3] < -2);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4]) >= 3;
