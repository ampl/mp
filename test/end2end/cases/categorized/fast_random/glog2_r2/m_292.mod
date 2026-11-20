var x{1..5} integer >= 0 <= 2;

maximize o: 5*x[1] - 3*x[2] + x[3] - 2*x[4] + 6*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 3;
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: (x[5] < 2) or (x[3] >= 1) or (x[2] >= 0);
