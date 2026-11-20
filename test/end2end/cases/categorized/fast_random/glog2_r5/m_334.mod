var x{1..3} integer >= -3 <= 1;
minimize o: (-3)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: exists {i in 1..3} x[i] == 0;
s.t. c3: (x[3] < 1) and (x[3] == -1) ==> (x[2] >= 1);
