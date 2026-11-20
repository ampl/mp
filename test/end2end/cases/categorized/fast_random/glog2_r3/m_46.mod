var x{1..4} integer >= 0 <= 4;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (if x[1] <= 1 then x[4] else x[2]) < 0;
s.t. c3: (if x[1] <= 0 then x[2] else x[4]) < 0;
s.t. c4: exists {i in 1..4} x[i] == 2;
