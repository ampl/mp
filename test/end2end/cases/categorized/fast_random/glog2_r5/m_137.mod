var x{1..6} integer >= 0 <= 3;
minimize o: (-3000)*x[1] + (-333)*x[2] + (51)*x[3] + (51)*x[4] + (-333)*x[5] + (-3000)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: abs(x[6] - x[2]) <= 3;
s.t. c3: (x[3] >= 3) ==> (x[1] < 2) else (x[1] != 1);
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
