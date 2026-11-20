var x{1..5} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (5)*x[2] + (4)*x[3] + (2)*x[4] + (0)*x[5];
s.t. c1: alldiff(x[4], x[5], x[3]);
s.t. c2: atleast 2 {i in 1..5} (x[i] >= 1);
s.t. c3: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) != 0;
s.t. c4: x[4] == 0;
s.t. c5: x[1] == 0;
