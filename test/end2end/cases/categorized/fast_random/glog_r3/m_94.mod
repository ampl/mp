var x{1..5} integer >= -2 <= 2;
minimize o: (3)*x[1] + (4)*x[2] + (1)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: (if x[3] < 0 then x[2] else x[1]) <= 2;
s.t. c2: abs(x[4] - x[1]) > 0;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
