var x{1..7} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (-5)*x[2] + (0)*x[3] + (-2)*x[4] + (1)*x[5] + (-3)*x[6] + (-4)*x[7];
s.t. c1: (x[7] < 1) and (x[5] <= 1) ==> (x[7] >= 1);
s.t. c2: forall {i in 1..6} x[i] != x[i+1];
s.t. c3: alldiff(x[6], x[1], x[2]);
