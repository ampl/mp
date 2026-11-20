var x{1..5} integer >= 0 <= 4;
minimize o: (4)*x[1] + (-1)*x[2] + (1)*x[3] + (3)*x[4] + (0)*x[5];
s.t. c1: ((if x[1] != 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c2: (x[3] != 5) and (x[3] != 5) ==> (x[5] >= 5);
s.t. c3: (x[3] > 5) ==> (x[1] != 0);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
