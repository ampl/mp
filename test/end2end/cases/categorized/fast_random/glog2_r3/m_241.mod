var x{1..6} integer >= 0 <= 4;
maximize o: (1)*x[1] + (5)*x[2] + (2)*x[3] + (-4)*x[4] + (1)*x[5] + (4)*x[6];
s.t. c1: abs(x[1] - x[1]) < 1;
s.t. c2: (if x[5] < 4 then x[3] else x[5]) != -1;
s.t. c3: alldiff(x[6], x[3]);
