var x{1..3} integer >= 0 <= 3;
maximize o: (5)*x[1] + (2)*x[2] + (2)*x[3];
s.t. c1: alldiff(x[1], x[3]);
s.t. c2: max(x[2], x[1]) - min(x[3], x[2]) < 1;
s.t. c3: (x[3] != 3) <==> (x[2] + x[2] >= 1);
s.t. c4: (if x[2] >= 2 then x[2] else x[3]) <= 2;
