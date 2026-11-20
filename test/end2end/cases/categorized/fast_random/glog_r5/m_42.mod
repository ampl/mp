var x{1..4} integer >= 0 <= 4;
maximize o: (3)*x[1] + (0)*x[2] + (-5)*x[3] + (-2)*x[4];
s.t. c1: (x[3] < 4) and (x[4] < 4) ==> (x[4] >= 4);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (if x[2] > 4 then x[3] else x[2]) < 2;
s.t. c4: ((if x[1] < 4 then 1 else 0) + (if x[3] >= 4 then 1 else 0)) >= 1;
