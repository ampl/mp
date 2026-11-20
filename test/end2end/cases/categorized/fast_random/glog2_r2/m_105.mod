var x{1..5} integer >= -2 <= 2;
maximize o: (-2)*x[1] + (6)*x[2] + (-1)*x[3] + (6)*x[4] + (4)*x[5];
s.t. c1: (x[5] != 2) ==> (x[3] <= -2) else (x[3] > 0);
s.t. c2: (if x[4] <= 1 then x[5] else x[2]) <= 2;
s.t. c3: atmost 3 {i in 1..5} (x[i] >= -2);
s.t. c4: (x[2] <= 1) ==> (x[2] < -2) else (x[1] >= 0);
s.t. c5: exists {i in 1..5} x[i] == 2;
s.t. c6: (x[2] == 0) <==> (x[1] + x[4] >= -1);
s.t. c7: not ((x[3] <= 1) and (x[1] <= 0));
