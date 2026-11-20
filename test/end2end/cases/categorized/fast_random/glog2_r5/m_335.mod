var x{1..4} integer >= 0 <= 2;

minimize o: -3*x[1] + 4*x[2] + 3*x[3] - x[4];

s.t. c1: max(x[2], x[1]) - min(x[3], x[3]) > 0;
s.t. c2: (if x[3] < 2 then x[3] else x[1]) != 2;
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: (x[2] > 0) ==> (x[2] > 1) else (x[2] <= 0);
s.t. c5: max(x[3], x[4]) - min(x[4], x[3]) >= 0;
s.t. c6: (x[1] < 0) ==> (x[3] >= 1);
