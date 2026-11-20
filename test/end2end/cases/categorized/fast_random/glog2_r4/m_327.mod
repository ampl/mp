var x{1..4} integer >= -2 <= 1;

maximize o: 51*x[1] - 666*x[2] - 333*x[3] + 3000*x[4];

s.t. c1: atmost 2 {i in 1..4} (x[i] >= -2);
s.t. c2: (if x[3] == 0 then x[2] else x[2]) <= 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) == 1;
s.t. c4: (if x[4] < -1 then x[2] else x[1]) != -1;
