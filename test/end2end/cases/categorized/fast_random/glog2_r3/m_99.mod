var x{1..6} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (2)*x[2] + (-5)*x[3] + (4)*x[4] + (1)*x[5] + (-4)*x[6];
s.t. c1: not ((x[3] != 3) and (x[5] != 2));
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: alldiff(x[1], x[3]);
s.t. c4: not ((x[3] != 3) and (x[5] != 2));
s.t. c5: alldiff(x[1], x[3]);
