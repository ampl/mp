var x{1..3} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (4)*x[2] + (-1)*x[3];
s.t. c1: abs(x[1] - x[1]) <= 0;
s.t. c2: x[3] == 1;
