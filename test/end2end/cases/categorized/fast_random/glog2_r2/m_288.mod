var x{1..6} integer >= 0 <= 1;
minimize o: (-2)*x[1] + (-4)*x[2] + (5)*x[3] + (6)*x[4] + (3)*x[5] + (1)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[2] <= 1) and (x[6] > 0) ==> (x[4] >= 1);
s.t. c3: not ((x[1] < 1) and (x[2] == 0));
