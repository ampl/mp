var x{1..5} integer >= 0 <= 3;

maximize o: -4*x[1] + x[2] + 0*x[3] - 5*x[4] + 3*x[5];

s.t. c1: max(x[3], x[1]) - min(x[2], x[3]) > 2;
s.t. c2: x[3] == 0;
