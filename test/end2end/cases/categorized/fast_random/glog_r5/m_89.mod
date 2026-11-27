var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] + 2*x[2] + 2*x[3];

s.t. c1: not ((x[2] == 2) and (x[3] == 1));
s.t. c2: (x[2] > 2) or (x[3] != 3) or (x[2] >= 3);
