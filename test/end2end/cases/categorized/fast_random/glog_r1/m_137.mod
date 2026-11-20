var x{1..6} integer >= -2 <= 0;
minimize o: (-2)*x[1] + (-5)*x[2] + (4)*x[3] + (4)*x[4] + (3)*x[5] + (2)*x[6];
s.t. c1: alldiff(x[2], x[5], x[4], x[1]);
s.t. c2: abs(x[1] - x[3]) >= 1;
