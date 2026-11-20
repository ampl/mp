var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] + 3*x[2] + 5*x[3] + 6*x[4] + 6*x[5] + x[6];

s.t. c1: not ((x[2] >= 3) and (x[6] < 0));
s.t. c2: max(x[6], x[3]) - min(x[5], x[2]) > 3;
s.t. c3: not ((x[2] >= 3) and (x[6] < 0));
s.t. c4: max(x[6], x[3]) - min(x[5], x[2]) > 3;
