var x{1..3} integer >= -3 <= 1;
maximize o: (-2)*x[1] + (-3)*x[2] + (2)*x[3];
s.t. c1: not ((x[1] >= -2) and (x[2] < 1));
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 0);
