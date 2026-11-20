var x{1..7} integer >= 0 <= 1;
minimize o: (-1)*x[1] + (-1)*x[2] + (-2)*x[3] + (-4)*x[4] + (-1)*x[5] + (4)*x[6] + (-1)*x[7];
s.t. c1: forall {i in 1..6} x[i] != x[i+1];
s.t. c2: ((if x[4] <= 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
