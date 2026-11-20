var x{1..4} integer >= 0 <= 2;
minimize o: (2)*x[1] + (3)*x[2] + (-4)*x[3] + (-5)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) == 3;
s.t. c2: max(x[2], x[4]) - min(x[3], x[1]) == 0;
s.t. c3: (x[1] >= 0) <==> (x[4] + x[1] != 2);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: forall {i in 1..3} x[i] != x[i+1];
