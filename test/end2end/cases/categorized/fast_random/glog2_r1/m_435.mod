var x{1..4} integer >= -2 <= 1;

minimize o: 6*x[1] + 2*x[2] + 5*x[3] + x[4];

s.t. c1: (x[1] <= 0) or (x[2] >= 0) or (x[4] < 0);
s.t. c2: not ((x[3] >= 0) and (x[4] <= -2));
s.t. c3: atmost 2 {i in 1..4} (x[i] >= 1);
s.t. c4: numberof -2 in (x[1], x[2], x[3], x[4]) != 1;
s.t. c5: (x[3] != 1) and (x[2] > -1) ==> (x[2] >= 1);
s.t. c6: (x[3] != -1) ==> (x[2] == 0) else (x[2] < -2);
