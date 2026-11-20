var x{1..5} integer >= 0 <= 3;
minimize o: (2)*x[1] + (-5)*x[2] + (4)*x[3] + (2)*x[4] + (3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[3] != 3) or (x[2] >= 3) or (x[1] > 2);
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: ((if x[1] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
