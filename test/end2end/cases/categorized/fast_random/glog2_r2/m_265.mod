var x{1..6} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (5)*x[2] + (1)*x[3] + (4)*x[4] + (-2)*x[5] + (5)*x[6];
s.t. c1: atmost 2 {i in 1..6} (x[i] >= 3);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 2;
s.t. c3: exists {i in 1..6} x[i] == 2;
s.t. c4: abs(x[4] - x[2]) >= 3;
s.t. c5: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 3;
