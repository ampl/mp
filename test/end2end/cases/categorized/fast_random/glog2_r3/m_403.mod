var x{1..4} integer >= -2 <= 1;
minimize o: (-5)*x[1] + (6)*x[2] + (-5)*x[3] + (-1)*x[4];
s.t. c1: abs(x[3] - x[1]) >= 2;
s.t. c2: alldiff(x[1], x[2], x[3], x[4]);
s.t. c3: (if x[3] < 0 then x[2] else x[1]) <= -1;
s.t. c4: x[1] == 1;
s.t. c5: x[3] == 0;
