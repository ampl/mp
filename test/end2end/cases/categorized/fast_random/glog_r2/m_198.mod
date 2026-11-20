var x{1..4} integer >= 0 <= 2;
maximize o: (-3)*x[1] + (0)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: not ((x[1] <= 0) and (x[4] <= 2));
s.t. c4: ((x[2] != 2) ==> (x[1] >= 1)) and ((x[1] >= 1) ==> (x[3] >= 2));
