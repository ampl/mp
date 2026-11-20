var x{1..6} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-3)*x[2] + (-1)*x[3] + (3)*x[4] + (1)*x[5] + (0)*x[6];
s.t. c1: atleast 2 {i in 1..6} (x[i] >= 2);
s.t. c2: x[6] == 1;
