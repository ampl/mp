var x{1..4} integer >= -1 <= 1;
minimize o: (4)*x[1] + (1)*x[2] + (5)*x[3] + (5)*x[4];
s.t. c1: (abs(x[2] - x[1]) != 2) or (abs(x[2] - x[1]) != 2);
s.t. c2: (if x[3] == x[2] then x[3] else x[2]) >= 1;
s.t. c3: (x[2] + x[4] > -1) <==> (x[3] == x[2]);
