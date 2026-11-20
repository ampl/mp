var x{1..5} integer >= -2 <= 1;
maximize o: (0)*x[1] + (-5)*x[2] + (5)*x[3] + (-2)*x[4] + (0)*x[5];
s.t. c1: max(x[3], x[3]) - min(x[3], x[5]) <= -2;
s.t. c2: (if x[5] <= 1 then x[4] else x[1]) <= 1;
s.t. c3: abs(x[5] - x[3]) < 2;
s.t. c4: (x[2] > 0) and (x[1] == -1) ==> (x[5] >= 0);
s.t. c5: x[5] == -2;
