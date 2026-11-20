var x{1..3} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (3)*x[2] + (-1)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: ((if x[2] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
