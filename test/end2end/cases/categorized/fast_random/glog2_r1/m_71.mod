var x{1..3} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (4)*x[2] + (-5)*x[3];
s.t. c1: not ((x[3] > 4) and (x[1] != 5));
s.t. c2: exists {i in 1..3} x[i] == 5;
s.t. c3: (x[3] <= 5) or (x[3] == -1) or (x[1] < 0);
