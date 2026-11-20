var x{1..6} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (-4)*x[2] + (3)*x[3] + (-4)*x[4] + (-1)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 3;
