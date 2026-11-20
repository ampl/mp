var x{1..4} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (5)*x[2] + (5)*x[3] + (4)*x[4];
s.t. c1: (if x[1] == 0 then x[2] else x[2]) != 4;
s.t. c2: (x[1] > 1) ==> (x[1] >= 4) else (x[3] != 2);
