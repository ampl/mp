var x{1..5} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: not ((x[1] >= 0) and (x[5] >= 2));
