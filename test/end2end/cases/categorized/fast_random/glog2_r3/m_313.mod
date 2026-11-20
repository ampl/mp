var x{1..3} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-1)*x[2] + (-5)*x[3];
s.t. c1: not ((x[2] != 1) and (x[2] >= 1));
s.t. c2: (if x[3] > 0 then x[3] else x[2]) < 0;
s.t. c3: abs(x[3] - x[1]) == 1;
s.t. c4: not ((x[2] < 0) and (x[3] < 0));
s.t. c5: x[2] == 2;
s.t. c6: x[3] == 1;
