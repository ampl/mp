var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (6)*x[2] + (-4)*x[3] + (3)*x[4] + (-3)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[2], x[5]) == 1;
s.t. c2: ((if x[3] < 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: (x[2] <= 0) <==> (x[2] + x[2] == 3);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
