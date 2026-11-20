var x{1..5} integer >= -2 <= 1;
minimize o: (-2)*x[1] + (1)*x[2] + (-3)*x[3] + (0)*x[4] + (-4)*x[5];
s.t. c1: max(x[4], x[2]) - min(x[2], x[4]) >= 1;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c3: (x[4] > 1) and (x[1] <= -2) ==> (x[4] >= 1);
s.t. c4: not ((x[5] >= -2) and (x[3] < 0));
s.t. c5: alldiff(x[4], x[2]);
s.t. c6: abs(x[4] - x[5]) < 1;
