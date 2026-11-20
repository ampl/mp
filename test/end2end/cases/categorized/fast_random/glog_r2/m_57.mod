var x{1..6} integer >= -2 <= 2;
maximize o: (-2)*x[1] + (-1)*x[2] + (-3)*x[3] + (-2)*x[4] + (-5)*x[5] + (-4)*x[6];
s.t. c1: not ((x[5] >= -1) and (x[2] <= -1));
s.t. c2: (x[5] < -1) ==> (x[6] == 2) else (x[5] >= 0);
s.t. c3: max(x[1], x[1]) - min(x[1], x[5]) > 0;
s.t. c4: exists {i in 1..6} x[i] == 2;
