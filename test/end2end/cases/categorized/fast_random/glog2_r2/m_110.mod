var x{1..4} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (-4)*x[2] + (2)*x[3] + (2)*x[4];
s.t. c1: (x[3] < 1) <==> (x[3] + x[1] <= 0);
s.t. c2: (if x[3] != 0 then x[1] else x[1]) < 2;
