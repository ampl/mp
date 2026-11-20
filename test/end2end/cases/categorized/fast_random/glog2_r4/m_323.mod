var x{1..3} integer >= -2 <= 0;
maximize o: (3)*x[1] + (5)*x[2] + (-4)*x[3];
s.t. c1: ((if x[3] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: x[2] == -1;
s.t. c3: x[1] == 0;
