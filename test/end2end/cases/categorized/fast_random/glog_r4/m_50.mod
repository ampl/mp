var x{1..4} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (1)*x[2] + (2)*x[3] + (-5)*x[4];
s.t. c1: not ((x[3] > 0) and (x[1] == 1));
s.t. c2: (if x[2] >= 1 then x[4] else x[2]) != 0;
s.t. c3: (if x[3] < 1 then x[3] else x[3]) <= 1;
s.t. c4: not ((x[1] == 1) and (x[2] > 1));
