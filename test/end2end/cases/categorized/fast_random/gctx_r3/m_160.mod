var x{1..4} integer >= -1 <= 1;
minimize o: (-3)*x[1] + (1)*x[2] + (2)*x[3] + (6)*x[4];
s.t. c1: (x[1] != 1) <==> (min(x[2], x[2]) != 1);
s.t. c2: not (x[1] != 1);
