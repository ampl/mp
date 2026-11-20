var x{1..5} integer >= -2 <= 1;
maximize o: (5)*x[1] + (1)*x[2] + (3)*x[3] + (0)*x[4] + (-2)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: ((if x[1] <= 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[1] < 0) ==> (x[1] == -2) else (x[4] != -1);
s.t. c4: (x[2] <= -1) ==> (x[1] <= -2);
s.t. c5: not ((x[5] == 1) and (x[4] <= -1));
s.t. c6: (x[2] == 1) and (x[1] <= -1) ==> (x[1] >= 1);
s.t. c7: abs(x[3] - x[1]) == 0;
