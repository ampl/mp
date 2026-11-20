var x{1..4} integer >= -1 <= 1;
minimize o: (5)*x[1] + (2)*x[2] + (-2)*x[3] + (1)*x[4] + (-6)*(if x[3] != 0 then 1 else 0);
s.t. c1: (max(x[4], x[2]) != 0) <==> (x[1] >= -1);
s.t. c2: (if x[1] >= -1 then x[3] else x[4]) == 1;
