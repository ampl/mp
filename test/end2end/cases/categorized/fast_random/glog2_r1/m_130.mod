var x{1..4} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (5)*x[2] + (0)*x[3] + (-3)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: max(x[1], x[3]) - min(x[4], x[4]) <= -2;
s.t. c3: abs(x[2] - x[1]) < 3;
