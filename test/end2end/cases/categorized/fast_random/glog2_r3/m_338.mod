var x{1..7} integer >= 0 <= 1;
minimize o: (2)*x[1] + (2)*x[2] + (2)*x[3] + (-4)*x[4] + (4)*x[5] + (-5)*x[6] + (4)*x[7];
s.t. c1: max(x[7], x[7]) - min(x[3], x[7]) >= 0;
s.t. c2: (if x[5] < 1 then x[1] else x[1]) == 0;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6], x[7]) != 3;
s.t. c4: abs(x[3] - x[7]) >= 0;
