var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (-5)*x[2] + (6)*x[3];
s.t. c1: not ((x[3] >= 0) and (x[1] == 0));
s.t. c2: (x[2] < 0) and (x[2] >= 1) ==> (x[2] >= 0);
s.t. c3: x[3] == 2;
