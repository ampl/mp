var x{1..3} integer >= 0 <= 3;
minimize o: (6)*x[1] + (1)*x[2] + (-3)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: ((x[2] > 3) ==> (x[3] < 3)) and ((x[3] < 3) ==> (x[3] >= 3));
s.t. c3: abs(x[2] - x[3]) > 1;
