var x{1..3} integer >= -2 <= 1;
minimize o: (3)*x[1] + (-3)*x[2] + (4)*x[3];
s.t. c1: not ((x[2] <= -3) and (x[1] > -2));
s.t. c2: (x[3] != 1) or (x[2] < -2) or (x[3] <= 1);
s.t. c3: not ((x[1] == 2) and (x[2] >= 2));
