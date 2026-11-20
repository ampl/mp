var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (5)*x[2] + (6)*x[3] + (6)*x[4];
s.t. c1: atmost 1 {i in 1..4} (x[i] >= 3);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: abs(x[4] - x[2]) != 0;
