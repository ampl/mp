var x{1..4} integer >= 0 <= 4;
minimize o: (6)*x[1] + (-5)*x[2] + (3)*x[3] + (2)*x[4];
s.t. c1: ((if x[2] > 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: x[3] == 2;
