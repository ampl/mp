var x{1..5} integer >= -3 <= 1;
maximize o: (-4)*x[1] + (-2)*x[2] + (-1)*x[3] + (1)*x[4] + (6)*x[5];
s.t. c1: abs(x[2] - x[5]) >= 3;
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: (x[1] >= -2) or (x[5] > -3) or (x[2] < -2);
