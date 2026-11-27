var x{1..6} integer >= -2 <= 0;

minimize o: 0*x[1] - 3*x[2] - 4*x[3] + 4*x[4] - 4*x[5] - x[6];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 3;
s.t. c3: exists {i in 1..6} x[i] == -2;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
s.t. c5: not ((x[5] <= 0) and (x[2] != -1));
s.t. c6: abs(x[1] - x[6]) == 0;
