var x{1..6} integer >= -2 <= 0;
maximize o: (4)*x[1] + (0)*x[2] + (2)*x[3] + (-2)*x[4] + (-2)*x[5] + (-4)*x[6];
s.t. c1: abs(x[4] - x[6]) < 0;
s.t. c2: ((if x[5] == 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: ((if x[1] != -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c4: x[5] == -1;
