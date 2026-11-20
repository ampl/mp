var x{1..6} integer >= 0 <= 4;
maximize o: (2)*x[1] + (0)*x[2] + (-5)*x[3] + (2)*x[4] + (1)*x[5] + (3)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: ((if x[6] >= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c3: abs(x[3] - x[4]) == 2;
s.t. c4: (x[6] < 4) or (x[1] > 0) or (x[1] == 3);
s.t. c5: atleast 3 {i in 1..6} (x[i] >= 1);
s.t. c6: alldiff(x[3], x[1], x[6], x[2]);
