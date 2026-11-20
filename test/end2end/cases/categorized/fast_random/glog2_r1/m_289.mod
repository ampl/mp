var x{1..3} integer >= 0 <= 2;
maximize o: (1)*x[1] + (1)*x[2] + (-4)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: max(x[3], x[1]) - min(x[3], x[2]) > 2;
s.t. c3: abs(x[2] - x[2]) >= 1;
s.t. c4: exists {i in 1..3} x[i] == 2;
s.t. c5: abs(x[2] - x[2]) >= 1;
