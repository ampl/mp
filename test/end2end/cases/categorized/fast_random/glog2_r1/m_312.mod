var x{1..4} integer >= 0 <= 2;
minimize o: (4)*x[1] + (5)*x[2] + (-4)*x[3] + (1)*x[4];
s.t. c1: abs(x[1] - x[1]) < 1;
s.t. c2: atleast 3 {i in 1..4} (x[i] >= 1);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c4: atleast 2 {i in 1..4} (x[i] >= 0);
