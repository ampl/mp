var x{1..3} integer >= -2 <= 1;

minimize o: 6*x[1] + 0*x[2] - 5*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= -2);
s.t. c2: ((if x[3] > -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c3: not ((x[3] >= -1) and (x[1] == -1));
s.t. c4: numberof -1 in (x[1], x[2], x[3]) == 3;
