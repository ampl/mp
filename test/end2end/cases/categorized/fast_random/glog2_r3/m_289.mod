var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (0)*x[2] + (5)*x[3] + (6)*x[4] + (0)*x[5];
s.t. c1: abs(x[2] - x[3]) > 2;
s.t. c2: (x[3] == 0) <==> (x[2] + x[2] == 0);
s.t. c3: not ((x[1] >= 1) and (x[3] >= 2));
s.t. c4: (if x[1] < 0 then x[2] else x[2]) != 1;
s.t. c5: abs(x[2] - x[3]) > 2;
s.t. c6: (if x[1] < 0 then x[2] else x[2]) != 1;
