var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] + 5*x[2] + 4*x[3];

s.t. c1: not ((x[2] < 4) and (x[2] < 4));
s.t. c2: not ((x[3] >= 3) and (x[1] <= 2));
s.t. c3: not ((x[2] < 4) and (x[2] < 4));
s.t. c4: not ((x[3] >= 3) and (x[1] <= 2));
