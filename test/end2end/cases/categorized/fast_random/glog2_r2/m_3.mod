var x{1..5} integer >= -2 <= 1;
maximize o: (1)*x[1] + (5)*x[2] + (6)*x[3] + (2)*x[4] + (6)*x[5];
s.t. c1: (if x[3] < 1 then x[5] else x[5]) != -2;
s.t. c2: abs(x[1] - x[5]) <= 2;
s.t. c3: not ((x[2] <= 1) and (x[4] > -2));
s.t. c4: x[2] == 1;
s.t. c5: x[5] == -2;
