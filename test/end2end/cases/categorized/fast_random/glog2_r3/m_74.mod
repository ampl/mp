var x{1..3} integer >= -3 <= 2;

maximize o: -2*x[1] - 4*x[2] + 3*x[3];

s.t. c1: abs(x[1] - x[1]) <= 5;
