var x{1..6} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-1)*x[2] + (4)*x[3] + (4)*x[4] + (5)*x[5] + (-1)*x[6];
s.t. c1: (x[4] > 0) ==> (x[5] >= 0);
s.t. c2: (if x[4] == 1 then x[1] else x[3]) != 1;
s.t. c3: max(x[3], x[3]) - min(x[5], x[5]) != 2;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
