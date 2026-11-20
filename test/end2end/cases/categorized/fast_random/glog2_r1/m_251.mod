var x{1..6} integer >= -2 <= 2;

maximize o: 3*x[1] + 6*x[2] - x[3] + 2*x[4] + 4*x[5] - 5*x[6];

s.t. c1: not ((x[2] >= -2) and (x[1] > -2));
