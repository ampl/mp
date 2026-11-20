var x{1..3} integer >= 0 <= 2;
minimize o: (-666)*x[1] + (-1000)*x[2] + (17)*x[3];
s.t. c1: (if x[1] == 1 then x[3] else x[1]) >= 1;
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
