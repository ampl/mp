var x{1..5} integer >= -2 <= 2;
maximize o: (0)*x[1] + (-2)*x[2] + (1)*x[3] + (3)*x[4] + (2)*x[5];
s.t. c1: atleast 3 {i in 1..5} (x[i] >= -2);
