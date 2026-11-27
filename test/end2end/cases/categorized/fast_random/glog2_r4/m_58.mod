var x{1..4} integer >= 0 <= 1;

maximize o: -3*x[1] + 4*x[2] - 4*x[3] - 5*x[4];

s.t. c1: abs(x[3] - x[4]) >= 0;
