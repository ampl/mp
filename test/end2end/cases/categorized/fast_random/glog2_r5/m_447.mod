var x{1..3} integer >= -2 <= 1;

minimize o: -3000*x[1] - 333*x[2] - 999*x[3];

s.t. c1: max(x[2], x[1]) - min(x[1], x[1]) != 0;
s.t. c2: ((x[2] <= 0) ==> (x[3] >= -2)) and ((x[3] >= -2) ==> (x[1] >= 0));
s.t. c3: not ((x[2] != -1) and (x[1] > 1));
s.t. c4: exists {i in 1..3} x[i] == -1;
