var x{1..3} integer >= -2 <= 1;

maximize o: 4*x[1] - x[2] + 3*x[3];

s.t. c1: not ((x[2] <= -1) and (x[3] != -1));
