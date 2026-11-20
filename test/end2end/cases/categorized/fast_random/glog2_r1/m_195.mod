var x{1..5} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (-1)*x[2] + (5)*x[3] + (-3)*x[4] + (6)*x[5];
s.t. c1: ((if x[3] > 4 then 1 else 0) + (if x[5] >= 4 then 1 else 0)) >= 1;
s.t. c2: (x[5] > 2) and (x[4] < 1) ==> (x[5] >= 2);
s.t. c3: exists {i in 1..5} x[i] == 2;
s.t. c4: abs(x[5] - x[4]) > 1;
s.t. c5: max(x[2], x[2]) - min(x[5], x[3]) >= 1;
s.t. c6: (x[3] == 4) and (x[1] != 3) ==> (x[4] >= 4);
