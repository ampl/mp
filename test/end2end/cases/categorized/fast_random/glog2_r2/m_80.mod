var x{1..4} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (0)*x[2] + (-1)*x[3] + (-2)*x[4];
s.t. c1: (x[1] > 1) and (x[1] != 0) ==> (x[2] >= 1);
s.t. c2: (if x[3] == 0 then x[2] else x[1]) >= -1;
