var x{1..5} integer >= -2 <= 2;

maximize o: 17*x[1] - 999*x[2] + 3000*x[3] + 1000*x[4] + 1000*x[5];

s.t. c1: abs(x[1] - x[3]) <= 4;
