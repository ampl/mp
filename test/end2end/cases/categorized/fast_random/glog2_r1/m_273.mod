var x{1..3} integer >= -2 <= 0;
maximize o: (-3)*x[1] + (-5)*x[2] + (-1)*x[3];
s.t. c1: max(x[1], x[3]) - min(x[2], x[3]) >= 0;
s.t. c2: x[1] == 0;
s.t. c3: x[3] == -2;
