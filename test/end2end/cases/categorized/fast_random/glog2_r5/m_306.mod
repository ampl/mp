var x{1..3} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: max(x[1], x[2]) - min(x[1], x[1]) < 0;
s.t. c2: ((if x[1] == 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
