var x{1..5} integer >= -2 <= 1;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: abs(x[2] - x[5]) <= 3;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: atleast 1 {i in 1..5} (x[i] >= 1);
