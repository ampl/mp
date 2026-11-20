var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] - 2*x[2] - 2*x[3] + 6*x[4];

s.t. c1: abs(x[2] - x[2]) == 3;
