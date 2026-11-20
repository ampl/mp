var x{1..3} integer >= -2 <= 0;
minimize o: (5)*x[1] + (6)*x[2] + (-4)*x[3];
s.t. c1: (x[2] != -2) ==> (x[2] != -1);
s.t. c2: (x[1] == -2) or (x[1] <= -1) or (x[2] > 0);
s.t. c3: abs(x[2] - x[3]) <= 0;
s.t. c4: exists {i in 1..3} x[i] == -2;
