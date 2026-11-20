var x{1..4} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (-4)*x[2] + (0)*x[3] + (-4)*x[4];
s.t. c1: max(x[3], x[2]) - min(x[2], x[1]) != 4;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: exists {i in 1..4} x[i] == 4;
s.t. c4: numberof 4 in (x[1], x[2], x[3], x[4]) == 1;
