var x{1..4} integer >= -2 <= 0;
maximize o: (-5)*x[1] + (1)*x[2] + (6)*x[3] + (-4)*x[4];
s.t. c1: not ((x[4] < -1) and (x[2] == 0));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[3] <= -1) and (x[4] > 0) ==> (x[1] >= -1);
s.t. c4: (if x[3] <= -2 then x[1] else x[4]) == 0;
