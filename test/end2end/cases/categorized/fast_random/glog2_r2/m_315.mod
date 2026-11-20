var x{1..3} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-1)*x[2] + (-4)*x[3];
s.t. c1: (x[2] != 0) ==> (x[3] < 1) else (x[3] >= 1);
s.t. c2: (x[1] >= 0) <==> (x[2] + x[3] <= 2);
s.t. c3: alldiff(x[2], x[1], x[3]);
s.t. c4: forall {i in 1..2} x[i] != x[i+1];
s.t. c5: forall {i in 1..2} x[i] != x[i+1];
