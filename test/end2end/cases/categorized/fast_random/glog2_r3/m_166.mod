var x{1..3} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: abs(x[2] - x[1]) == 2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: numberof 1 in (x[1], x[2], x[3]) == 0;
s.t. c4: ((x[3] != 1) ==> (x[3] < 1)) and ((x[3] < 1) ==> (x[3] >= 1));
