var x{1..6} integer >= 0 <= 4;
minimize o: (-1)*x[1] + (5)*x[2] + (0)*x[3] + (5)*x[4] + (5)*x[5] + (4)*x[6];
s.t. c1: not ((x[4] < 1) and (x[6] < 3));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
