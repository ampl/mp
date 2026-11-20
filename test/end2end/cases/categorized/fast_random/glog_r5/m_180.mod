var x{1..3} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (4)*x[2] + (2)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[1] >= 0) or (x[3] == 3) or (x[3] >= 3);
