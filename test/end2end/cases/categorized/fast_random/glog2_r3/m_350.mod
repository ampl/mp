var x{1..4} integer >= -2 <= 1;
maximize o: (-3)*x[1] + (1)*x[2] + (1)*x[3] + (-4)*x[4];
s.t. c1: ((if x[1] > 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: abs(x[4] - x[2]) <= 2;
