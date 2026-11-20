var x{1..3} integer >= 0 <= 2;
maximize o: (3)*x[1] + (2)*x[2] + (-2)*x[3];
s.t. c1: atleast 1 {i in 1..3} (x[i] >= 1);
