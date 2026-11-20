var x{1..5} integer >= -2 <= 0;

maximize o: 5*x[1] - 3*x[2] - 5*x[3] + 3*x[4] + 6*x[5];

s.t. c1: abs(x[3] - x[1]) != 0;
s.t. c2: abs(x[3] - x[1]) != 0;
