var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: abs(x[5] - x[3]) == 2;
s.t. c2: (x[4] > 2) <==> (x[5] + x[4] > 3);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
