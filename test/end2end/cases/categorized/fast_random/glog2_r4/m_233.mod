var x{1..6} integer >= 0 <= 2;
minimize o: (3)*x[1] + (6)*x[2] + (-3)*x[3] + (-3)*x[4] + (2)*x[5] + (5)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= 1);
s.t. c2: (x[5] <= 2) ==> (x[6] > 1);
s.t. c3: forall {i in 1..5} x[i] != x[i+1];
s.t. c4: x[2] == 1;
