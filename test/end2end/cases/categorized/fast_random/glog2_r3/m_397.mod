var x{1..5} integer >= 0 <= 3;

maximize o: 51*x[1] + 34*x[2] + 3000*x[3] + 51*x[4] + 51*x[5];

s.t. c1: not ((x[4] > 0) and (x[4] != 1));
