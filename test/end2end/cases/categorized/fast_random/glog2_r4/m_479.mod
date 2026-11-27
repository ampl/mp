var x{1..3} integer >= -2 <= 1;

maximize o: 4*x[1] - x[2] - 4*x[3];

s.t. c1: (x[3] == 1) or (x[2] == -1) or (x[3] <= -2);
s.t. c2: (x[3] == 1) or (x[2] == -1) or (x[3] <= -2);
