var x{1..5} integer >= -3 <= 2;
minimize o: (-5)*x[1] + (3)*x[2] + (-3)*x[3] + (0)*x[4] + (2)*x[5];
s.t. c1: atmost 3 {i in 1..5} (x[i] >= -3);
