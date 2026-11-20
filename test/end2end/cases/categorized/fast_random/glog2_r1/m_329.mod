var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (1)*x[2] + (4)*x[3] + (-4)*x[4] + (5)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[1] > 0) or (x[2] > 2) or (x[3] < 0);
s.t. c3: abs(x[3] - x[1]) >= 0;
s.t. c4: exists {i in 1..5} x[i] == 1;
s.t. c5: abs(x[3] - x[1]) >= 0;
