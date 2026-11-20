var x{1..3} integer >= -2 <= 1;

maximize o: 3*x[1] - 2*x[2] - 5*x[3];

s.t. c1: abs(x[1] - x[3]) == 2;
s.t. c2: abs(x[1] - x[3]) == 2;
