var x{1..4} integer >= 0 <= 2;
maximize o: (-1000)*x[1] + (3000)*x[2] + (-333)*x[3] + (3000)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[1] > 0) ==> (x[2] > 1) else (x[4] != 0);
s.t. c3: not ((x[3] <= 1) and (x[2] == 1));
