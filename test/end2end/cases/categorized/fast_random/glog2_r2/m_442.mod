var x{1..3} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-5)*x[2] + (2)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: ((if x[2] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
