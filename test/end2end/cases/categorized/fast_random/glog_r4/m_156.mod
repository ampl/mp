var x{1..6} integer >= -2 <= 1;

maximize o: 5*x[1] + 6*x[2] - 5*x[3] + 2*x[4] - 3*x[5] - 5*x[6];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 0;
s.t. c2: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c3: (if x[5] != -1 then x[4] else x[2]) >= -2;
s.t. c4: atmost 3 {i in 1..6} (x[i] >= 1);
