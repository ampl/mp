var x{1..3} integer >= 0 <= 4;
maximize o: (5)*x[1] + (3)*x[2] + (5)*x[3];
s.t. c1: alldiff(x[2], x[1], x[3]);
s.t. c2: (x[1] > 4) ==> (x[3] < 0) else (x[2] != 2);
s.t. c3: (x[3] != 3) ==> (x[2] != 2);
s.t. c4: (if x[1] >= 2 then x[1] else x[3]) > 0;
