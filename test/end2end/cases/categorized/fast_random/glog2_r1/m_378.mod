var x{1..7} integer >= 0 <= 1;
minimize o: (-1)*x[1] + (6)*x[2] + (-3)*x[3] + (-4)*x[4] + (2)*x[5] + (0)*x[6] + (-5)*x[7];
s.t. c1: forall {i in 1..6} x[i] != x[i+1];
