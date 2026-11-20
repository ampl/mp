var x{1..4} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (5)*x[2] + (0)*x[3] + (-3)*x[4];
s.t. c1: atmost 2 {i in 1..4} (x[i] >= 0);
s.t. c2: (x[2] >= 0) <==> (x[3] + x[1] <= 2);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: ((if x[2] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
