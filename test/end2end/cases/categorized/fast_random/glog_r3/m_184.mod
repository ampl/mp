var x{1..5} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (1)*x[2] + (-4)*x[3] + (-1)*x[4] + (1)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[2] > 3) ==> (x[4] >= 2);
