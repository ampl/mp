var x{1..5} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-4)*x[2] + (-4)*x[3] + (-2)*x[4] + (5)*x[5];
s.t. c1: atmost 3 {i in 1..5} (x[i] >= 1);
