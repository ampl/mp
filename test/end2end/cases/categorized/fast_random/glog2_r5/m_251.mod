var x{1..4} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (3)*x[2] + (1)*x[3] + (4)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
