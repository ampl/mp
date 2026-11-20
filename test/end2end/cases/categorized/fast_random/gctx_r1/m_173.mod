var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (6)*x[2] + (5)*x[3];
s.t. c1: (if x[1] + x[1] < 1 then 1 else 0) + (if x[3] + x[1] > 2 then 1 else 0) + (if x[1] + x[1] < 1 then 1 else 0) >= 1;
s.t. c2: (x[3] + x[1] > 2) ==> ((x[2] >= x[1]) or (x[2] >= x[1]));
