var x{1..4} integer >= 0 <= 4;
maximize o: (-5)*x[1] + (-4)*x[2] + (-3)*x[3] + (-5)*x[4];
s.t. c1: not ((x[2] >= 1) and (x[4] <= 3));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: atmost 1 {i in 1..4} (x[i] >= 4);
s.t. c4: (x[4] == 0) ==> (x[3] > 1) else (x[1] <= 0);
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4]) >= 1;
s.t. c6: max(x[2], x[4]) - min(x[3], x[1]) > 3;
s.t. c7: (x[1] <= 4) or (x[4] == 4) or (x[3] == 2);
