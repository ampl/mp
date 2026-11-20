var x{1..6} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-5)*x[2] + (3)*x[3] + (3)*x[4] + (-2)*x[5] + (2)*x[6];
s.t. c1: max(x[4], x[1]) - min(x[1], x[4]) == 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
s.t. c3: abs(x[5] - x[2]) >= 0;
s.t. c4: (x[4] <= 2) ==> (x[6] >= 2);
s.t. c5: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
