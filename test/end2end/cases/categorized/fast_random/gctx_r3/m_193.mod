var x{1..3} integer >= 0 <= 3;
minimize o: (6)*x[1] + (5)*x[2] + (0)*x[3];
s.t. c1: (x[1] == 0) or (min(x[3], x[2]) > 0);
s.t. c2: (max(x[3], x[2]) != 2) <==> (min(x[3], x[2]) > 0);
s.t. c3: not ((min(x[3], x[2]) > 0) and (x[1] == 0));
