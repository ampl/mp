var x{1..3} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (2)*x[2] + (4)*x[3];
s.t. c1: (x[2] >= 4) <==> (x[1] + x[1] >= 2);
s.t. c2: abs(x[1] - x[2]) >= 3;
s.t. c3: not ((x[3] < 3) and (x[1] != 2));
