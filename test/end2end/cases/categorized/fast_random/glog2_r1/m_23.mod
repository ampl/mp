var x{1..6} integer >= 0 <= 4;
maximize o: (4)*x[1] + (6)*x[2] + (-2)*x[3] + (4)*x[4] + (5)*x[5] + (1)*x[6];
s.t. c1: ((if x[2] < 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
s.t. c2: x[2] == 3;
s.t. c3: x[5] == 2;
