var x{1..6} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (-1)*x[2] + (-5)*x[3] + (1)*x[4] + (3)*x[5] + (5)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
