var x{1..3} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-4)*x[2] + (5)*x[3];
s.t. c1: (x[2] >= 0) <==> (x[2] < x[3]);
s.t. c2: not ((x[2] < x[3]) and (x[2] >= 0));
