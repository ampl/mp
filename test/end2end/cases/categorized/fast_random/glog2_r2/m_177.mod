var x{1..3} integer >= 0 <= 3;
maximize o: (34)*x[1] + (1000)*x[2] + (17)*x[3];
s.t. c1: (x[1] == 3) or (x[3] == 1) or (x[2] >= 3);
s.t. c2: exists {i in 1..3} x[i] == 2;
