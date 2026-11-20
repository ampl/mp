var x{1..4} integer >= -2 <= 2;

maximize o: x[1] - 2*x[2] + 6*x[3] + 5*x[4];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4]) > 0;
s.t. c2: (x[4] == -1) or (x[4] < 1) or (x[1] == 1);
s.t. c3: exists {i in 1..4} x[i] == 2;
s.t. c4: (x[3] <= -1) ==> (x[2] >= -1) else (x[2] < 0);
s.t. c5: (x[3] >= 1) ==> (x[1] != 2);
s.t. c6: (x[1] > 1) and (x[1] > 0) ==> (x[4] >= 1);
s.t. c7: max(x[3], x[3]) - min(x[2], x[4]) > -1;
