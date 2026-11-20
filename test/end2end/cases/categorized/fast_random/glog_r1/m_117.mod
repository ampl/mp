var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] + 6*x[2] + 3*x[3] - x[4];

s.t. c1: (x[1] == 2) ==> (x[3] < 1) else (x[1] <= 1);
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: max(x[3], x[4]) - min(x[3], x[4]) == 4;
