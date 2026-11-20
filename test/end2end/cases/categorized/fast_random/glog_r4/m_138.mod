var x{1..3} integer >= 0 <= 2;
minimize o: (1)*x[1] + (2)*x[2] + (1)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: numberof 2 in (x[1], x[2], x[3]) != 3;
