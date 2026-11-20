var x{1..4} integer >= -1 <= 1;
maximize o: (2)*x[1] + (3)*x[2] + (-1)*x[3] + (-5)*x[4];
s.t. c1: (x[3] != x[2]) <==> (x[2] == 0);
s.t. c2: (x[3] != x[2]) or (x[2] <= -1);
s.t. c3: (x[3] != x[2]) <==> (x[3] != x[2]);
s.t. c4: (if x[3] != x[2] then x[1] else 0) + (if x[2] <= -1 then x[2] else 0) > 1;
