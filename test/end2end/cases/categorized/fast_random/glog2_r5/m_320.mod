var x{1..5} integer >= -2 <= 1;
minimize o: (0)*x[1] + (4)*x[2] + (4)*x[3] + (2)*x[4] + (6)*x[5];
s.t. c1: atleast 3 {i in 1..5} (x[i] >= 1);
s.t. c2: alldiff(x[4], x[3], x[5]);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
