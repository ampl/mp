var x{1..4} integer >= -2 <= 1;

minimize o: -3*x[1] - 3*x[2] + 6*x[3] - 3*x[4];

s.t. c1: max(x[3], x[4]) - min(x[2], x[4]) > -1;
s.t. c2: atmost 1 {i in 1..4} (x[i] >= -2);
s.t. c3: (x[1] == -2) <==> (x[1] + x[4] == -1);
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) >= 0;
s.t. c5: (x[4] >= -1) or (x[4] > 1) or (x[4] < 0);
s.t. c6: (x[3] == -2) ==> (x[3] > 0) else (x[2] < 1);
