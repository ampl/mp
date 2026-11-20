var x{1..4} integer >= -1 <= 1;
minimize o: (4)*x[1] + (-4)*x[2] + (-2)*x[3] + (3)*x[4] + (-6)*(if x[2] + x[1] == -1 then 1 else 0);
s.t. c1: (x[2] + x[1] == -1) <==> (x[2] + x[1] == -1);
s.t. c2: not ((abs(x[3] - x[4]) != 2) and (abs(x[3] - x[4]) != 2));
