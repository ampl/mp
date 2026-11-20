var x{1..5} integer >= 0 <= 1;
minimize o: (-4)*x[1] + (-2)*x[2] + (2)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: (x[2] >= 1) or (x[2] >= 0) or (x[3] <= 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[5] != 1) or (x[1] < 1) or (x[1] != 0);
s.t. c4: abs(x[5] - x[5]) > 1;
