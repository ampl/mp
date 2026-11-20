var x{1..3} integer >= 0 <= 2;
maximize o: (2)*x[1] + (4)*x[2] + (-1)*x[3];
s.t. c1: (if max(x[3], x[3]) < 0 then x[2] else x[1]) == 2;
s.t. c2: (max(x[3], x[3]) < 0) or (x[2] < 0);
s.t. c3: not ((x[3] < x[3]) and (max(x[3], x[3]) < 0));
s.t. c4: not ((max(x[3], x[3]) < 0) and (x[3] < x[3]));
