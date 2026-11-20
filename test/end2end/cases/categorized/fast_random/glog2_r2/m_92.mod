var x{1..3} integer >= -2 <= 2;
maximize o: (-5)*x[1] + (0)*x[2] + (4)*x[3];
s.t. c1: (x[3] < 2) and (x[3] < 1) ==> (x[3] >= 2);
