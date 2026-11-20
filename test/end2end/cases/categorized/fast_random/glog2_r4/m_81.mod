var x{1..3} integer >= 0 <= 4;
minimize o: (6)*x[1] + (6)*x[2] + (2)*x[3];
s.t. c1: (x[1] < 4) or (x[2] <= -1) or (x[3] == 5);
s.t. c2: abs(x[1] - x[1]) < 4;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
