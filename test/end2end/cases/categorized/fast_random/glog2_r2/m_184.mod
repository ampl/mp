var x{1..3} integer >= -3 <= 1;
minimize o: (5)*x[1] + (-1)*x[2] + (1)*x[3];
s.t. c1: (if x[2] < -2 then x[2] else x[1]) > -1;
s.t. c2: max(x[1], x[1]) - min(x[3], x[3]) < 1;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
