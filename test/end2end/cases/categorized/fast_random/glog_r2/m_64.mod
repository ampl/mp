var x{1..3} integer >= -2 <= 0;
minimize o: (-2)*x[1] + (1)*x[2] + (2)*x[3];
s.t. c1: numberof -2 in (x[1], x[2], x[3]) == 3;
s.t. c2: atleast 1 {i in 1..3} (x[i] >= -1);
