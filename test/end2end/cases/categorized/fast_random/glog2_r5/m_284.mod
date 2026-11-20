var x{1..5} integer >= -3 <= 1;
minimize o: (-1)*x[1] + (-3)*x[2] + (4)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) < 1;
s.t. c2: abs(x[2] - x[4]) >= 4;
s.t. c3: (x[4] <= -3) and (x[4] > -1) ==> (x[3] >= -3);
s.t. c4: exists {i in 1..5} x[i] == -2;
