var x{1..3} integer >= -3 <= 1;

maximize o: -3*x[1] + 3*x[2] + 0*x[3];

s.t. c1: abs(x[2] - x[2]) == 2;
