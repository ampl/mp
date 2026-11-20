var x{1..4} integer >= 0 <= 2;

maximize o: x[1] - x[2] - 4*x[3] - 5*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) > 1;
s.t. c2: (x[1] != 1) ==> (x[3] != 1);
s.t. c3: max(x[4], x[2]) - min(x[3], x[3]) <= 0;
s.t. c4: exists {i in 1..4} x[i] == 2;
s.t. c5: (x[3] > 1) or (x[2] > 0) or (x[4] <= 1);
s.t. c6: not ((x[4] <= 2) and (x[2] < 1));
