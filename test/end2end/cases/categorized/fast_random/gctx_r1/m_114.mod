var x{1..4} integer >= 0 <= 3;
minimize o: (1)*x[1] + (-3)*x[2] + (2)*x[3] + (5)*x[4];
s.t. c1: not ((abs(x[4] - x[2]) <= 2) and (x[4] > 2));
s.t. c2: (abs(x[4] - x[2]) <= 2) or (abs(x[3] - x[3]) <= 0);
s.t. c3: (x[4] > 2) and ((abs(x[4] - x[2]) <= 2) or not (abs(x[4] - x[2]) <= 2));
s.t. c4: (abs(x[4] - x[2]) <= 2) <==> (x[4] > 2);
