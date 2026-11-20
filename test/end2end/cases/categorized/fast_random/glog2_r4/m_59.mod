var x{1..4} integer >= 0 <= 4;
maximize o: (5)*x[1] + (-5)*x[2] + (2)*x[3] + (5)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[1] != 3) <==> (x[2] + x[1] > 1);
s.t. c3: (x[2] > 1) ==> (x[4] != 4) else (x[2] == 3);
s.t. c4: (x[4] != 4) or (x[2] != 0) or (x[2] != 2);
s.t. c5: forall {i in 1..3} x[i] != x[i+1];
s.t. c6: (x[4] != 4) or (x[2] != 0) or (x[2] != 2);
