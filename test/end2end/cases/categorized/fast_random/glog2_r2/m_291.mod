var x{1..3} integer >= -2 <= 2;
minimize o: (-2)*x[1] + (-3)*x[2] + (3)*x[3];
s.t. c1: numberof 3 in (x[1], x[2], x[3]) >= 2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: numberof -3 in (x[1], x[2], x[3]) != 1;
s.t. c4: ((if x[3] > -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
