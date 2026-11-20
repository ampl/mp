var x{1..3} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: max(x[3], x[3]) - min(x[1], x[1]) < 1;
s.t. c2: ((if x[1] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[3] <= 3) and (x[3] >= 4) ==> (x[1] >= 3);
