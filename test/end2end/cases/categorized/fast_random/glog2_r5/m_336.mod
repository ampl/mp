var x{1..5} integer >= -2 <= 1;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: (x[4] < 1) or (x[5] < -2) or (x[5] != 1);
s.t. c2: abs(x[5] - x[2]) != 1;
s.t. c3: exists {i in 1..5} x[i] == 1;
