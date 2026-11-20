var x{1..3} integer >= 0 <= 2;
maximize o: (3)*x[1] + (3)*x[2] + (-4)*x[3];
s.t. c1: (x[3] >= 0) <==> (x[3] + x[3] < 0);
s.t. c2: (if x[3] != 1 then x[3] else x[2]) > 1;
s.t. c3: not ((x[1] == 1) and (x[2] != 1));
