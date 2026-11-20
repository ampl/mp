var x{1..6} integer >= -2 <= 1;
minimize o: (3)*x[1] + (6)*x[2] + (2)*x[3] + (-4)*x[4] + (1)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: abs(x[4] - x[1]) == 0;
s.t. c3: abs(x[6] - x[5]) <= 0;
s.t. c4: exists {i in 1..6} x[i] == 1;
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: exists {i in 1..6} x[i] == 1;
