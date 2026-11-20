var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (-1)*x[2] + (3)*x[3];
s.t. c1: atleast 1 {i in 1..3} (x[i] >= 0);
