var x{1..4} integer >= 0 <= 4;

minimize o: 5*x[1] + 4*x[2] + 5*x[3] + 4*x[4];

s.t. c1: not ((x[4] < 3) and (x[2] > 4));
s.t. c2: not ((x[1] > 3) and (x[2] > 1));
s.t. c3: max(x[3], x[3]) - min(x[2], x[2]) < 1;
s.t. c4: not ((x[4] < 3) and (x[2] > 4));
s.t. c5: max(x[3], x[3]) - min(x[2], x[2]) < 1;
