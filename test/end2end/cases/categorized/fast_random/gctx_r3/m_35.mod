var x{1..4} integer >= -1 <= 1;
minimize o: (2)*x[1] + (2)*x[2] + (-5)*x[3] + (2)*x[4] + (-5)*(if max(x[3], x[2]) >= 0 then 1 else 0);
s.t. c1: (min(x[2], x[1]) != 0) <==> (x[4] != x[4]);
s.t. c2: (if min(x[2], x[1]) != 0 then x[2] else 0) + (if max(x[3], x[2]) >= 0 then x[4] else 0) >= 0;
s.t. c3: (max(x[1], x[1]) <= -1) <==> (x[4] != x[4]);
s.t. c4: not ((x[4] != x[4]) and (min(x[2], x[1]) != 0));
