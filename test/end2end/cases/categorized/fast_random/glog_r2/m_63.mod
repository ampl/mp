var x{1..6} integer >= 0 <= 2;
maximize o: (-5)*x[1] + (0)*x[2] + (-5)*x[3] + (-5)*x[4] + (-4)*x[5] + (2)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: not ((x[3] > 1) and (x[6] >= 1));
s.t. c3: alldiff(x[1], x[4], x[5], x[2]);
s.t. c4: (x[6] == 1) or (x[3] < 2) or (x[5] == 1);
