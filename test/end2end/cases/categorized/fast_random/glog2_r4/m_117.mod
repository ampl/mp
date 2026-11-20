var x{1..6} integer >= -2 <= 2;

maximize o: -1000*x[1] + 17*x[2] + 34*x[3] + 17*x[4] + 51*x[5] + 51*x[6];

s.t. c1: abs(x[1] - x[1]) != 0;
