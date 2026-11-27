var x{1..4} integer >= 0 <= 3;

maximize o: 4*x[1] + 3*x[2] + 6*x[3] - 3*x[4];

s.t. c1: not ((x[4] < 0) and (x[4] == 0));
s.t. c2: max(x[1], x[2]) - min(x[2], x[4]) != 1;
s.t. c3: max(x[4], x[3]) - min(x[3], x[4]) > 1;
s.t. c4: not ((x[4] < 0) and (x[4] == 0));
s.t. c5: max(x[4], x[3]) - min(x[3], x[4]) > 1;
