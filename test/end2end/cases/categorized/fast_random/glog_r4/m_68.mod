var x{1..5} integer >= -2 <= 1;
minimize o: (-1)*x[1] + (0)*x[2] + (1)*x[3] + (2)*x[4] + (-2)*x[5];
s.t. c1: atleast 3 {i in 1..5} (x[i] >= 1);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) != 1;
