var x{1..3} integer >= -3 <= 1;
minimize o: (-1)*x[1] + (-2)*x[2] + (-2)*x[3];
s.t. c1: (x[3] != -1) and (x[2] > -2) ==> (x[2] >= -1);
s.t. c2: max(x[2], x[2]) - min(x[2], x[1]) <= 1;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
