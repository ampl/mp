var x{1..4} integer >= -2 <= 0;

minimize o: -5*x[1] + 3*x[2] + 6*x[3] + 5*x[4];

s.t. c1: max(x[4], x[3]) - min(x[2], x[3]) >= -1;
s.t. c2: exists {i in 1..4} x[i] == -1;
s.t. c3: (x[2] <= -2) ==> (x[3] == -1) else (x[2] == -2);
s.t. c4: numberof -2 in (x[1], x[2], x[3], x[4]) == 1;
