var x{1..5} integer >= -2 <= 2;

maximize o: 3*x[1] - x[2] + 2*x[3] + 5*x[4] + x[5];

s.t. c1: (x[2] != 2) ==> (x[3] >= 1);
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 1;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) > 3;
