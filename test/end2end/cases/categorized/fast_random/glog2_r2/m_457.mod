var x{1..3} integer >= -2 <= 1;

maximize o: -3000*x[1] - 333*x[2] + 2000*x[3];

s.t. c1: not ((x[2] > -2) and (x[1] >= 1));
