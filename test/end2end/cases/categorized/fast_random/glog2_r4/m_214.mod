var x{1..5} integer >= -3 <= 2;

maximize o: -2*x[1] + x[2] - 4*x[3] - 2*x[4] + 2*x[5];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) > 3;
s.t. c2: atmost 3 {i in 1..5} (x[i] >= 0);
