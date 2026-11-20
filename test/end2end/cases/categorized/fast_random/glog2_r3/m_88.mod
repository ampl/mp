var x{1..6} integer >= 0 <= 1;
minimize o: (5)*x[1] + (-2)*x[2] + (5)*x[3] + (5)*x[4] + (2)*x[5] + (-1)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: max(x[6], x[3]) - min(x[2], x[2]) != 0;
s.t. c3: (x[1] >= 1) <==> (x[6] + x[1] == 0);
