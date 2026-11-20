var x{1..3} integer >= -1 <= 2;
maximize o: (2)*x[1] + (2)*x[2] + (2)*x[3];
s.t. c1: (x[3] >= 1) ==> (max(x[1], x[1]) <= -1);
s.t. c2: (x[1] > x[1]) ==> ((x[3] >= 1) or (min(x[3], x[3]) > 2));
s.t. c3: (min(x[3], x[3]) > 2) and ((min(x[3], x[3]) > 2) or not (max(x[1], x[1]) <= -1));
