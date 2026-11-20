var x{1..6} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (5)*x[2] + (4)*x[3] + (2)*x[4] + (1)*x[5] + (-4)*x[6];
s.t. c1: (x[1] >= 2) or (x[1] <= 1) or (x[1] > 2);
s.t. c2: (x[6] <= 2) <==> (x[1] + x[4] > 0);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 2;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
