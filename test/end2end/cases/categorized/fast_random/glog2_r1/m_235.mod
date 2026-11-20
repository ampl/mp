var x{1..5} integer >= 0 <= 2;
maximize o: (4)*x[1] + (-5)*x[2] + (6)*x[3] + (0)*x[4] + (-2)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) > 1;
s.t. c3: (x[4] == 0) or (x[5] > 1) or (x[1] <= 1);
s.t. c4: (x[5] == 0) ==> (x[2] > 0);
s.t. c5: alldiff(x[2], x[4], x[3]);
s.t. c6: (x[2] != 0) ==> (x[2] >= 0);
s.t. c7: (x[2] < 2) ==> (x[4] != 0) else (x[2] == 1);
