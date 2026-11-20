var x{1..6} integer >= 0 <= 2;

maximize o: x[1] + 4*x[2] + 6*x[3] + 4*x[4] + 2*x[5] + 2*x[6];

s.t. c1: atmost 2 {i in 1..6} (x[i] >= 0);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
s.t. c3: atmost 2 {i in 1..6} (x[i] >= 0);
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
