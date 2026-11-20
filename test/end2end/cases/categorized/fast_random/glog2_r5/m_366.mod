var x{1..5} integer >= -2 <= 0;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) != 1;
s.t. c3: ((if x[2] < -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
