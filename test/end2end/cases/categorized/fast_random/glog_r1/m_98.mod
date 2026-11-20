var x{1..5} integer >= 0 <= 3;

maximize o: 4*x[1] - x[2] - 5*x[3] - 3*x[4] + 6*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) >= 0;
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: (x[1] >= 2) <==> (x[3] + x[1] > 3);
