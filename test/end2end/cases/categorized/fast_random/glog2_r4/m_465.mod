var x{1..5} integer >= 0 <= 2;

minimize o: -2*x[1] + x[2] + x[3] + 5*x[4] - 4*x[5];

s.t. c1: (x[4] <= 2) ==> (x[1] == 0);
s.t. c2: max(x[4], x[2]) - min(x[2], x[5]) == 0;
s.t. c3: (x[1] >= 0) ==> (x[3] > 1) else (x[4] >= 0);
s.t. c4: exists {i in 1..5} x[i] == 2;
s.t. c5: abs(x[4] - x[4]) == 2;
s.t. c6: exists {i in 1..5} x[i] == 0;
