var x{1..4} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (-3)*x[2] + (2)*x[3] + (5)*x[4];
s.t. c1: (max(x[3], x[2]) != 2) <==> (x[3] >= x[4]);
s.t. c2: (max(x[3], x[2]) != 2) or (x[3] >= x[4]);
s.t. c3: not ((max(x[3], x[2]) != 2) and (max(x[3], x[2]) != 2));
s.t. c4: not (x[3] >= x[4]);
