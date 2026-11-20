var x{1..5} integer >= 0 <= 4;
maximize o: (4)*x[1] + (2)*x[2] + (-3)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (x[1] < 0) <==> (x[1] + x[1] > 4);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: alldiff(x[3], x[4], x[2], x[5]);
