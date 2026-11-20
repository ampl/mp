var x{1..4} integer >= 0 <= 4;

maximize o: -x[1] - 3*x[2] - 3*x[3] - 2*x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= -1);
s.t. c2: numberof 4 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c3: not ((x[3] != 0) and (x[2] == 5));
