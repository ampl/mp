var x{1..5} integer >= 0 <= 2;
minimize o: (6)*x[1] + (2)*x[2] + (-2)*x[3] + (-2)*x[4] + (1)*x[5];
s.t. c1: ((x[4] >= 1) ==> (x[1] <= 1)) and ((x[1] <= 1) ==> (x[3] >= 1));
s.t. c2: abs(x[1] - x[2]) != 2;
s.t. c3: max(x[3], x[4]) - min(x[5], x[3]) != 2;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
