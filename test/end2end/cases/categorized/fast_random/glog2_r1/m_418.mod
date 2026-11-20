var x{1..5} integer >= 0 <= 1;
maximize o: (-5)*x[1] + (2)*x[2] + (2)*x[3] + (-4)*x[4] + (-4)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: ((if x[3] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..5} x[i] == 0;
