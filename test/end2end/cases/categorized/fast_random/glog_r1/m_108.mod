var x{1..4} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-4)*x[2] + (2)*x[3] + (0)*x[4];
s.t. c1: (x[1] < 0) <==> (x[1] + x[4] >= 1);
s.t. c2: (if x[1] > 2 then x[1] else x[4]) != 1;
s.t. c3: not ((x[4] > 0) and (x[4] == 1));
