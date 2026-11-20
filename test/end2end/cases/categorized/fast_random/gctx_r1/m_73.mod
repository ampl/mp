var x{1..3} integer >= 0 <= 3;
minimize o: (3)*x[1] + (-3)*x[2] + (6)*x[3];
s.t. c1: (if x[3] >= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0) + (if x[2] != 3 then 1 else 0) >= 2;
s.t. c2: (x[3] == x[1]) ==> (x[2] != 3) else (x[3] == x[1]);
