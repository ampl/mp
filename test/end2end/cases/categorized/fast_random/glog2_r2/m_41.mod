var x{1..4} integer >= 0 <= 2;

minimize o: 5*x[1] - x[2] - 4*x[3] + 4*x[4];

s.t. c1: max(x[2], x[4]) - min(x[1], x[2]) != 0;
s.t. c2: (x[4] == 2) and (x[2] <= 2) ==> (x[1] >= 2);
s.t. c3: exists {i in 1..4} x[i] == 3;
