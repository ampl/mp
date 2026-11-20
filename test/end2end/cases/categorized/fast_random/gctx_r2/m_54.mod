var x{1..4} integer >= -1 <= 2;
maximize o: (4)*x[1] + (-5)*x[2] + (-5)*x[3] + (-5)*x[4];
s.t. c1: (if x[3] == x[3] then x[1] else 0) + (if x[3] == x[3] then x[3] else 0) < 1;
s.t. c2: not ((max(x[4], x[3]) == 1) and (x[3] == x[3]));
s.t. c3: (x[3] == x[3]) <==> (max(x[4], x[3]) == 1);
