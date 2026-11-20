var x{1..5} integer >= -3 <= 1;
minimize o: (3)*x[1] + (5)*x[2] + (-4)*x[3] + (0)*x[4] + (-1)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
