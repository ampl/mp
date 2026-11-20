var x{1..5} integer >= 0 <= 4;

maximize o: -3*x[1] - 4*x[2] - 3*x[3] + 4*x[4] + 6*x[5];

s.t. c1: atmost 2 {i in 1..5} (x[i] >= 2);
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
s.t. c3: not ((x[2] == 1) and (x[5] > 3));
