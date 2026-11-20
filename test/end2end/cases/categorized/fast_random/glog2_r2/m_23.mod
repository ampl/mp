var x{1..3} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (-3)*x[2] + (-5)*x[3];
s.t. c1: (x[3] < 3) or (x[3] > 3) or (x[2] < 4);
s.t. c2: x[1] == 2;
