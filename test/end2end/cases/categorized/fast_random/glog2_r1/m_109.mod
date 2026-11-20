var x{1..4} integer >= 0 <= 3;

minimize o: 4*x[1] + 5*x[2] - 4*x[3] - 2*x[4];

s.t. c1: max(x[2], x[1]) - min(x[3], x[4]) != 1;
s.t. c2: abs(x[2] - x[2]) != 0;
s.t. c3: abs(x[4] - x[2]) > 1;
s.t. c4: max(x[2], x[1]) - min(x[3], x[4]) != 1;
s.t. c5: abs(x[4] - x[2]) > 1;
