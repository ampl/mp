var x{1..4} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (3)*x[2] + (3)*x[3] + (-3)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: max(x[4], x[3]) - min(x[4], x[2]) != -2;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4]) == 1;
s.t. c4: (x[2] > 0) and (x[4] >= 2) ==> (x[3] >= 0);
