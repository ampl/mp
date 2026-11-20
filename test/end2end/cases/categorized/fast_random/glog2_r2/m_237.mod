var x{1..3} integer >= -2 <= 0;
minimize o: (51)*x[1] + (-2000)*x[2] + (1000)*x[3];
s.t. c1: not ((x[2] <= -1) and (x[1] != 0));
s.t. c2: numberof 0 in (x[1], x[2], x[3]) >= 3;
s.t. c3: exists {i in 1..3} x[i] == 0;
