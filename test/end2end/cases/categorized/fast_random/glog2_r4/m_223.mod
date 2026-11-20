var x{1..4} integer >= -2 <= 2;
minimize o: (4)*x[1] + (5)*x[2] + (0)*x[3] + (-2)*x[4];
s.t. c1: (x[2] == 1) <==> (x[3] + x[2] >= 0);
s.t. c2: (if x[2] != -2 then x[2] else x[4]) != 0;
s.t. c3: x[4] == -2;
