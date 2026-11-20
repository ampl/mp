var x{1..4} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (3)*x[2] + (2)*x[3] + (6)*x[4];
s.t. c1: abs(x[4] - x[4]) == 0;
s.t. c2: abs(x[4] - x[2]) > 0;
s.t. c3: (x[3] >= 1) ==> (x[2] > 0) else (x[4] <= 2);
s.t. c4: (x[1] == 2) ==> (x[3] == 1);
s.t. c5: exists {i in 1..4} x[i] == 0;
s.t. c6: exists {i in 1..4} x[i] == 2;
s.t. c7: not ((x[4] <= 0) and (x[2] <= 0));
