var x{1..4} integer >= 0 <= 3;
minimize o: (2)*x[1] + (5)*x[2] + (-4)*x[3] + (-3)*x[4];
s.t. c1: (x[3] != 1) and (x[1] < 0) ==> (x[3] >= 1);
s.t. c2: exists {i in 1..4} x[i] == 3;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: x[1] == 0;
