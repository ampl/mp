var x{1..3} integer >= 0 <= 3;

maximize o: 3*x[1] + 4*x[2] + x[3];

s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: numberof 2 in (x[1], x[2], x[3]) < 3;
s.t. c3: (x[2] > 3) or (x[2] < 2) or (x[2] == 2);
