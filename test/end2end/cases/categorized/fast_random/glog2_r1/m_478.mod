var x{1..4} integer >= 0 <= 1;
minimize o: (5)*x[1] + (-1)*x[2] + (3)*x[3] + (-5)*x[4];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c2: not ((x[4] > 1) and (x[3] >= 0));
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
