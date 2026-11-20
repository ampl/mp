var x{1..5} integer >= 0 <= 4;
maximize o: (5)*x[1] + (3)*x[2] + (3)*x[3] + (-5)*x[4] + (4)*x[5];
s.t. c1: ((if x[2] < 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[5] != 2) <==> (x[1] + x[5] >= 4);
s.t. c3: (x[4] <= 0) or (x[3] > 1) or (x[3] < 0);
s.t. c4: (x[1] == 3) ==> (x[1] == 0);
s.t. c5: (x[2] >= 2) and (x[5] >= 0) ==> (x[2] >= 2);
s.t. c6: exists {i in 1..5} x[i] == 0;
