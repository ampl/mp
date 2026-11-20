var x{1..5} integer >= -3 <= 2;
minimize o: (1)*x[1] + (-5)*x[2] + (-5)*x[3] + (-1)*x[4] + (6)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) > 1;
