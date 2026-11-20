var x{1..3} integer >= -2 <= 2;
minimize o: (-3)*x[1] + (3)*x[2] + (6)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: not ((x[1] >= -1) and (x[1] >= 0));
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: not ((x[3] <= 0) and (x[1] <= 1));
s.t. c5: x[3] == -1;
s.t. c6: x[1] == -1;
