var x{1..5} integer >= -2 <= 1;

maximize o: x[1] + 6*x[2] - 4*x[3] - 4*x[4] + x[5];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) == 0;
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: max(x[2], x[5]) - min(x[1], x[2]) <= 1;
s.t. c4: atmost 2 {i in 1..5} (x[i] >= 0);
s.t. c5: (x[1] < 1) or (x[4] == 0) or (x[4] != 0);
