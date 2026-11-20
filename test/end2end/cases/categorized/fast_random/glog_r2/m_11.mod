var x{1..3} integer >= 0 <= 2;

maximize o: -2*x[1] - 4*x[2] + 0*x[3];

s.t. c1: abs(x[2] - x[1]) >= 0;
s.t. c2: not ((x[3] < 2) and (x[1] < 0));
