var x{1..4} integer >= -2 <= 2;

minimize o: 6*x[1] + 5*x[2] + 5*x[3] + 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c3: (x[3] >= 2) ==> (x[2] > -1);
