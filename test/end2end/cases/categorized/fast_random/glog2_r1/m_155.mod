var x{1..6} integer >= 0 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (-4)*x[3] + (5)*x[4] + (-1)*x[5] + (-2)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: exists {i in 1..6} x[i] == 0;
s.t. c3: ((x[1] < 0) ==> (x[2] != 1)) and ((x[2] != 1) ==> (x[6] >= 0));
s.t. c4: (x[5] >= 2) <==> (x[2] + x[2] == 0);
s.t. c5: ((if x[6] < 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
