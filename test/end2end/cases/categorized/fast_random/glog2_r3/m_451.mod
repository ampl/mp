var x{1..3} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (2)*x[2] + (-4)*x[3];
s.t. c1: (x[1] <= 3) ==> (x[1] == 0);
