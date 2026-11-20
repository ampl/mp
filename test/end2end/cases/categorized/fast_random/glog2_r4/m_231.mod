var x{1..3} integer >= -2 <= 2;
minimize o: (2)*x[1] + (-2)*x[2] + (6)*x[3];
s.t. c1: (x[1] >= 3) or (x[1] <= 2) or (x[3] != 2);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: ((if x[2] < -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
