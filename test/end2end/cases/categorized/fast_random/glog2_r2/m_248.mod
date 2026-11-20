var x{1..4} integer >= 0 <= 1;
minimize o: (-1)*x[1] + (1)*x[2] + (5)*x[3] + (-5)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
