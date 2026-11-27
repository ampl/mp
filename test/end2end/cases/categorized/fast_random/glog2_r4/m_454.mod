var x{1..4} integer >= -2 <= 3;

maximize o: x[1] + 2*x[2] - x[3] + x[4];

s.t. c1: max(x[2], x[2]) - min(x[3], x[4]) <= 0;
s.t. c2: (x[4] > -1) or (x[2] > -2) or (x[3] == -1);
s.t. c3: max(x[2], x[3]) - min(x[3], x[4]) == 2;
