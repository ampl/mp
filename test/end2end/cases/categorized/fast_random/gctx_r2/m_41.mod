var x{1..4} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-4)*x[2] + (-5)*x[3] + (6)*x[4];
s.t. c1: (max(x[1], x[3]) < 2) <==> (max(x[1], x[3]) < 2);
s.t. c2: (x[4] + x[2] != 0) and ((x[2] > 0) or not (x[4] + x[2] != 0));
s.t. c3: (max(x[1], x[3]) < 2) and ((x[4] == x[4]) or not (max(x[1], x[3]) < 2));
