var x{1..6} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (3)*x[2] + (-5)*x[3] + (0)*x[4] + (1)*x[5] + (1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: abs(x[1] - x[5]) == 0;
s.t. c3: not ((x[3] == 3) and (x[2] <= 4));
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 1;
s.t. c5: exists {i in 1..6} x[i] == 2;
s.t. c6: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 1;
