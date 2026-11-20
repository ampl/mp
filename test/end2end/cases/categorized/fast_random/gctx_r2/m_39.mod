var x{1..4} integer >= 0 <= 2;
minimize o: (0)*x[1] + (3)*x[2] + (5)*x[3] + (2)*x[4];
s.t. c1: (x[4] == 2) and ((x[4] != x[2]) or not (x[4] == 2));
s.t. c2: (x[4] == 2) ==> (x[4] != x[2]) else (x[4] != x[2]);
