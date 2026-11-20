var x{1..4} integer >= 0 <= 3;
minimize o: (5)*x[1] + (-3)*x[2] + (4)*x[3] + (4)*x[4] + (0)*(if x[1] == 3 then 1 else 0);
s.t. c1: (if x[2] != 2 then x[2] else 0) + (if x[2] != 2 then x[4] else 0) <= 1;
s.t. c2: (x[1] == 3) ==> ((x[4] > x[2]) or (x[2] != 2));
