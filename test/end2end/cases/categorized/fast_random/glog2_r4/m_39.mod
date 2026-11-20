var x{1..4} integer >= 0 <= 4;

minimize o: -5*x[1] - 5*x[2] - 4*x[3] + 6*x[4];

s.t. c1: (x[1] == 4) or (x[1] > 0) or (x[1] != 1);
s.t. c2: not ((x[4] < 2) and (x[1] < 0));
s.t. c3: (x[1] == 4) or (x[1] > 0) or (x[1] != 1);
s.t. c4: not ((x[4] < 2) and (x[1] < 0));
