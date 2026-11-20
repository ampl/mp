var x{1..3} integer >= -2 <= 3;
minimize o: (-1)*x[1] + (3)*x[2] + (6)*x[3];
s.t. c1: (x[2] == 1) ==> (x[1] == 3) else (x[2] > 1);
s.t. c2: ((if x[3] >= 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
