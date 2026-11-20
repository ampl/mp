var x{1..4} integer >= 0 <= 4;

maximize o: -5*x[1] + 3*x[2] + 5*x[3] - x[4];

s.t. c1: ((if x[1] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c3: max(x[4], x[1]) - min(x[4], x[2]) > 4;
s.t. c4: atmost 1 {i in 1..4} (x[i] >= 1);
s.t. c5: x[1] == 0;
s.t. c6: x[4] == 3;
