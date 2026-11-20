var x{1..5} integer >= -1 <= 1;
minimize o: (5)*x[1] + (-3)*x[2] + (1)*x[3] + (-1)*x[4] + (3)*x[5] + (5)*(if max(x[1], x[4]) == 1 then 1 else 0);
s.t. c1: (if max(x[1], x[4]) == 1 then 1 else 0) + (if x[2] > 1 then 1 else 0) + (if max(x[2], x[4]) <= -1 then 1 else 0) >= 1;
s.t. c2: (max(x[1], x[4]) == 1) ==> (max(x[1], x[4]) == 1) else (max(x[2], x[4]) <= -1);
