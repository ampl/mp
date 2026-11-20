var x{1..6} integer >= -2 <= 1;
maximize o: (5)*x[1] + (-3)*x[2] + (3)*x[3] + (0)*x[4] + (-3)*x[5] + (0)*x[6];
s.t. c1: ((if x[3] == -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c2: abs(x[6] - x[4]) > 0;
s.t. c3: ((if x[3] == -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c4: abs(x[6] - x[4]) > 0;
