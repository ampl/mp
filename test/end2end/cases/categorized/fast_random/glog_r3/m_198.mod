var x{1..4} integer >= 0 <= 4;
maximize o: (2)*x[1] + (2)*x[2] + (-1)*x[3] + (-2)*x[4];
s.t. c1: not ((x[1] == 4) and (x[1] == 4));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[4] > 0) <==> (x[4] + x[2] == 4);
