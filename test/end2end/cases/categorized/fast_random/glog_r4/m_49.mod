var x{1..5} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (2)*x[2] + (6)*x[3] + (3)*x[4] + (-1)*x[5];
s.t. c1: (x[4] > 1) or (x[1] >= 0) or (x[4] > 2);
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: (x[2] >= 0) or (x[4] <= 1) or (x[3] <= 0);
