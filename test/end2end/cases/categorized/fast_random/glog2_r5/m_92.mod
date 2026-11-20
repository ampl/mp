var x{1..6} integer >= 0 <= 4;
minimize o: (6)*x[1] + (-5)*x[2] + (-1)*x[3] + (-4)*x[4] + (0)*x[5] + (-2)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
