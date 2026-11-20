var x{1..6} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-1)*x[2] + (5)*x[3] + (6)*x[4] + (-1)*x[5] + (-1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: (if x[1] <= 1 then x[6] else x[6]) == 0;
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: abs(x[4] - x[6]) > 1;
s.t. c5: not ((x[4] >= 1) and (x[6] > 2));
s.t. c6: ((x[2] == 2) ==> (x[1] < 1)) and ((x[1] < 1) ==> (x[6] >= 2));
s.t. c7: (x[1] == 2) ==> (x[1] != 2) else (x[3] >= 0);
