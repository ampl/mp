var x{1..3} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (-2)*x[2] + (1)*x[3];
s.t. c1: (x[3] >= 0) ==> (x[3] != 4) else (x[3] < 0);
s.t. c2: abs(x[3] - x[3]) > 3;
s.t. c3: exists {i in 1..3} x[i] == 3;
s.t. c4: max(x[3], x[1]) - min(x[1], x[3]) < 4;
s.t. c5: (if x[2] != 2 then x[3] else x[2]) != 3;
s.t. c6: abs(x[2] - x[3]) == 0;
