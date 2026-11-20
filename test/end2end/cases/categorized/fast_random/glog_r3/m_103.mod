var x{1..6} integer >= 0 <= 2;

maximize o: 5*x[1] + x[2] + x[3] + 2*x[4] - 3*x[5] + 6*x[6];

s.t. c1: abs(x[5] - x[5]) <= 2;
s.t. c2: not ((x[6] > 2) and (x[2] < 2));
