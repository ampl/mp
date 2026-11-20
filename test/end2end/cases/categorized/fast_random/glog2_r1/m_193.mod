var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] + 2*x[2] - 3*x[3] - 2*x[4] + 3*x[5];

s.t. c1: not ((x[4] >= 3) and (x[5] != 1));
s.t. c2: (x[5] <= 0) ==> (x[1] <= 1) else (x[5] > 2);
s.t. c3: max(x[2], x[2]) - min(x[1], x[3]) <= 3;
s.t. c4: exists {i in 1..5} x[i] == 0;
s.t. c5: x[4] == 0;
s.t. c6: x[5] == 2;
