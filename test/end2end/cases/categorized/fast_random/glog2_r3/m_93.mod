var x{1..5} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-3)*x[2] + (5)*x[3] + (4)*x[4] + (1)*x[5];
s.t. c1: (x[5] >= 1) and (x[2] >= 1) ==> (x[5] >= 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: exists {i in 1..5} x[i] == 0;
s.t. c4: x[1] == 2;
s.t. c5: x[4] == 0;
