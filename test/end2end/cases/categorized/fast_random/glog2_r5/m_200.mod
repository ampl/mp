var x{1..5} integer >= 0 <= 3;
minimize o: (3)*x[1] + (0)*x[2] + (1)*x[3] + (6)*x[4] + (1)*x[5];
s.t. c1: ((x[2] > 0) ==> (x[4] != 3)) and ((x[4] != 3) ==> (x[3] >= 0));
s.t. c2: atleast 2 {i in 1..5} (x[i] >= 3);
s.t. c3: abs(x[5] - x[5]) != 0;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
