var x{1..4} integer >= -2 <= 2;

minimize o: -2*x[1] - 3*x[2] - 3*x[3] + 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 2;
s.t. c2: (x[2] == -1) and (x[4] < 2) ==> (x[3] >= -1);
s.t. c3: (x[2] > 2) ==> (x[4] == 1) else (x[1] <= 0);
s.t. c4: not ((x[3] == 1) and (x[4] > 1));
s.t. c5: exists {i in 1..4} x[i] == 2;
s.t. c6: not ((x[3] == 1) and (x[4] > 1));
