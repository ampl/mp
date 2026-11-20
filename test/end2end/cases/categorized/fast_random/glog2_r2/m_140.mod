var x{1..5} integer >= -2 <= 0;
minimize o: (5)*x[1] + (0)*x[2] + (-1)*x[3] + (-3)*x[4] + (2)*x[5];
s.t. c1: not ((x[5] != -2) and (x[4] >= -1));
s.t. c2: (x[3] < 0) ==> (x[4] >= -1) else (x[4] < -1);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
