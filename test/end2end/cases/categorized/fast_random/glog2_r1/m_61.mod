var x{1..5} integer >= 0 <= 4;

maximize o: 3*x[1] - 3*x[2] - 5*x[3] + x[4] - 4*x[5];

s.t. c1: (x[2] > 4) and (x[2] >= 0) ==> (x[1] >= 4);
s.t. c2: exists {i in 1..5} x[i] == 5;
s.t. c3: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) > 0;
