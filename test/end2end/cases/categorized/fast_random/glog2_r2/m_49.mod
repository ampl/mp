var x{1..3} integer >= -2 <= 0;
maximize o: (0)*x[1] + (0)*x[2] + (-1)*x[3];
s.t. c1: exists {i in 1..3} x[i] == -2;
s.t. c2: exists {i in 1..3} x[i] == -2;
