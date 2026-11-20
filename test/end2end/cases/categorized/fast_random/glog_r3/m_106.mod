var x{1..5} integer >= -2 <= 0;
minimize o: (4)*x[1] + (-2)*x[2] + (1)*x[3] + (3)*x[4] + (-4)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
s.t. c3: exists {i in 1..5} x[i] == -1;
