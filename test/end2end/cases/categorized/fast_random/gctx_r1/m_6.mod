var x{1..4} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-1)*x[2] + (5)*x[3] + (-5)*x[4];
s.t. c1: (x[2] == 1) or (x[2] == 1);
s.t. c2: (x[2] == 1) <==> (max(x[1], x[3]) < 3);
s.t. c3: (if max(x[1], x[3]) < 3 then x[1] else x[2]) == 3;
s.t. c4: not ((max(x[1], x[3]) < 3) and (x[2] == 1));
