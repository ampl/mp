var x{1..6} integer >= 0 <= 4;
maximize o: (0)*x[1] + (4)*x[2] + (2)*x[3] + (2)*x[4] + (0)*x[5] + (-5)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: (x[5] == 1) <==> (x[1] + x[6] >= 3);
s.t. c3: exists {i in 1..6} x[i] == 4;
s.t. c4: (x[6] < 1) ==> (x[2] < 1);
s.t. c5: (x[6] < 4) <==> (x[4] + x[1] <= 4);
s.t. c6: (if x[6] >= 0 then x[5] else x[2]) > 3;
