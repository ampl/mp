var x{1..3} integer >= 0 <= 3;
maximize o: (3)*x[1] + (4)*x[2] + (1)*x[3];
s.t. c1: not (x[1] < 0);
s.t. c2: (if x[1] < 0 then x[1] else x[3]) == 0;
s.t. c3: not (x[1] < 0);
