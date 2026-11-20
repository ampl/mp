var x{1..3} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (6)*x[2] + (-2)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: (x[3] > 4) ==> (x[1] < 3);
s.t. c3: abs(x[1] - x[1]) == 1;
