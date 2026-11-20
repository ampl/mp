var x{1..3} integer >= -2 <= 0;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: (x[1] != -2) or (x[3] < -2) or (x[2] == 0);
s.t. c2: ((x[1] > 0) ==> (x[1] >= 0)) and ((x[1] >= 0) ==> (x[2] >= 0));
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
