var x{1..5} integer >= -2 <= 3;
maximize o: (-1)*x[1] + (-2)*x[2] + (-1)*x[3] + (3)*x[4] + (3)*x[5];
s.t. c1: (x[4] == 3) or (x[1] == 2) or (x[5] > 3);
s.t. c2: (x[2] != 1) or (x[1] >= 3) or (x[3] == 2);
s.t. c3: (if x[1] >= 2 then x[2] else x[2]) <= -2;
