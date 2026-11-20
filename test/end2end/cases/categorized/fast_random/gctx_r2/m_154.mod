var x{1..5} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (2)*x[2] + (-2)*x[3] + (3)*x[4] + (4)*x[5] + (2)*(if max(x[2], x[2]) >= 1 then 1 else 0);
s.t. c1: not (max(x[2], x[2]) >= 1);
s.t. c2: (min(x[2], x[3]) != 2) <==> (max(x[2], x[2]) >= 1);
