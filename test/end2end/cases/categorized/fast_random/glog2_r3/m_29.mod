var x{1..6} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-1)*x[2] + (-5)*x[3] + (-2)*x[4] + (1)*x[5] + (0)*x[6];
s.t. c1: alldiff(x[5], x[3]);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: max(x[1], x[4]) - min(x[5], x[2]) == 3;
s.t. c4: alldiff(x[5], x[3]);
s.t. c5: max(x[1], x[4]) - min(x[5], x[2]) == 3;
