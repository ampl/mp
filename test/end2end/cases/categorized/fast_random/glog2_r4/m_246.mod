var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: ((if x[1] == 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3]) >= 0;
s.t. c3: exists {i in 1..3} x[i] == 2;
