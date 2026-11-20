var x{1..6} integer >= 0 <= 3;
minimize o: (0)*x[1] + (4)*x[2] + (4)*x[3] + (-5)*x[4] + (3)*x[5] + (0)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
