var x{1..3} integer >= -1 <= 2;
maximize o: (3)*x[1] + (4)*x[2] + (2)*x[3];
s.t. c1: not (x[2] != x[3]);
s.t. c2: (x[2] != x[3]) <==> (x[2] + x[3] >= 1);
s.t. c3: not ((x[2] != x[3]) and (x[1] > -1));
s.t. c4: not (x[2] + x[3] >= 1);
