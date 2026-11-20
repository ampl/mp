var x{1..5} integer >= 0 <= 3;
minimize o: (4)*x[1] + (-5)*x[2] + (-3)*x[3] + (4)*x[4] + (2)*x[5];
s.t. c1: not ((x[2] > 0) and (x[2] > 0));
s.t. c2: (x[2] > 0) or (x[2] > 0);
