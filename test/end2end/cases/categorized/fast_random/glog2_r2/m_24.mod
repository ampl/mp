var x{1..3} integer >= -3 <= 1;
minimize o: (-5)*x[1] + (2)*x[2] + (-1)*x[3];
s.t. c1: abs(x[2] - x[1]) <= 1;
s.t. c2: exists {i in 1..3} x[i] == 1;
