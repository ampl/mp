var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[4] < 0) ==> (x[4] != 2);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) > 1;
s.t. c3: exists {i in 1..4} x[i] == 2;
