var x{1..4} integer >= 0 <= 2;

maximize o: 3*x[1] + 4*x[2] + 2*x[3] - x[4];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) > 3;
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: max(x[4], x[3]) - min(x[4], x[1]) != 1;
s.t. c4: not ((x[4] <= 0) and (x[3] < 2));
s.t. c5: (x[2] <= 2) ==> (x[1] <= 2);
