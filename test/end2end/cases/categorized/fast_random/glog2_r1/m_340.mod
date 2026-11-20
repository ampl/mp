var x{1..5} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (0)*x[2] + (4)*x[3] + (6)*x[4] + (5)*x[5];
s.t. c1: abs(x[5] - x[5]) != 2;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[5] != 0) ==> (x[2] >= 2);
s.t. c4: max(x[4], x[4]) - min(x[2], x[1]) > 0;
